"""
Author: Marco Maneta
Email: mmaneta@ekiconsult.com
"""

import math
from functools import singledispatchmethod
from logging import getLogger

import fiona
import flopy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from flopy.mf6.mfbase import MFDataException
from flopy.utils import GridIntersect
from flopy.utils import Raster
from flopy.utils.lgrutil import Lgr
from rasterio import features
from shapely.geometry import MultiLineString, shape

logger = getLogger(__name__)

REACH_FLAG = 999999


class NestedDomain:
    def __init__(self,
                 sim=None,
                 gwf=None,
                 modelgrid=None):

        self.sim = sim
        self.gwf = gwf

        self.parent_model_name = None

        if gwf is not None:
            modelgrid = self.gwf.modelgrid
            self.parent_model_name = gwf.name
        self.parent_domain = modelgrid

        self.lst_subdomain_names = []
        self.lst_subdomain_lgr = []

    @classmethod
    def from_parent_model(cls, gwf_model_name: str = None, **kwargs):
        try:
            sim = flopy.mf6.MFSimulation.load(**kwargs)
            gwf = sim.get_model(gwf_model_name)
        except MFDataException as e:
            raise FileNotFoundError

        return cls(sim, gwf)

    def define_subdomain(self, name: str,
                         istart: int = None,
                         istop: int = None,
                         jstart: int = None,
                         jstop: int = None,
                         kstart: int = None,
                         kstop: int = None,
                         nested_domain_shp: str = None,
                         feature_field: str = 'id',
                         xoff: float = 0.0,
                         yoff: float = 0.0,
                         angrot: float = 0.0,
                         num_cells_per_parent_cell: int = 3,
                         num_layers_per_parent_layer: list = 1):

        if istart == istop == jstart == jstop == kstart == kstop == nested_domain_shp is None:
            print("Either the bounding box or a shapefile to define the subdomain is required")
            return

        self.lst_subdomain_names.append(name)
        idomain = self.gwf.dis.idomain.get_data()

        if nested_domain_shp is None:
            # deactivate the cells in the parent domain where the child grid will be placed
            idomain[kstart:kstop + 1, istart:istop + 1, jstart:jstop + 1] = 2
        else:
            pols = fiona.open(nested_domain_shp)
            self.gwf.modelgrid.set_coord_info(xoff=xoff, yoff=yoff, angrot=angrot, crs=pols.crs)
            rst_tpl = Raster.raster_from_array(np.zeros((self.gwf.dis.nrow.data,
                                                         self.gwf.dis.ncol.data)),
                                               crs=pols.crs,
                                               modelgrid=self.gwf.modelgrid)

            for f in pols:
                shapes = (f.geometry, f.properties[feature_field])
                burned = features.rasterize(shapes=shapes,
                                            fill=0,
                                            out_shape=(self.gwf.dis.nrow.data,
                                                       self.gwf.dis.ncol.data),
                                            transform=rst_tpl.transform,
                                            all_touched=True)
                ir, ic = burned.nonzero()
                idomain[kstart:kstop, ir, ic] = 2

        self.lst_subdomain_lgr.append(
            Lgr(nlayp=self.gwf.dis.nlay.data,
                nrowp=self.gwf.dis.nrow.data,
                ncolp=self.gwf.dis.ncol.data,
                delrp=self.gwf.dis.delr.data,
                delcp=self.gwf.dis.delc.data,
                topp=self.gwf.dis.top.data,
                botmp=self.gwf.dis.botm.data,
                idomainp=idomain,  # self.gwf.dis.idomain.data,
                ncpp=num_cells_per_parent_cell,
                ncppl=num_layers_per_parent_layer,
                xllp=self.gwf.modelgrid.xoffset,
                yllp=self.gwf.modelgrid.yoffset,
                )
        )

        self.gwf.dis.idomain.set_data(self.lst_subdomain_lgr[-1].parent.idomain)

        # self.gwf.dis.idomain.set_data(idomain)
        self._display_domain_info()

    def _create_exchange_data(self, lgr, subdomain_name: str):
        logger.info("Creating exchange data for subdomain {}. This can take a minute...".format(subdomain_name))
        exchangedata = lgr.get_exchange_data(angldegx=True, cdist=True)
        exg = flopy.mf6.ModflowGwfgwf(self.sim,
                                      exgtype="GWF6-GWF6",
                                      xt3d=True,
                                      auxiliary=["angldegx", "cdist"],
                                      exgmnamea=self.gwf.name,
                                      exgmnameb=subdomain_name,
                                      nexg=len(exchangedata),
                                      exchangedata=exchangedata,
                                      pname=f"{subdomain_name}.exg",
                                      )
        return exg

    def get_flow_simulation(self):

        lst_exchg = []
        for lgr, name in zip(self.lst_subdomain_lgr, self.lst_subdomain_names):
            lst_exchg.append(
                self._create_exchange_data(lgr, name)
            )
        return NestedDomainSimulation(self.sim,
                                      self.parent_model_name,
                                      self.lst_subdomain_names,
                                      self.lst_subdomain_lgr)

    def _display_domain_info(self):
        for i, d in enumerate(self.lst_subdomain_names):

            print(f"DOMAIN {d}:")
            for name, p in zip(["Parent grid", "Nested grid"],
                               [self.lst_subdomain_lgr[i].parent, self.lst_subdomain_lgr[i].child]):
                print(f"\t{d}: {name}")
                print(f"\t\tNum layers {p.nlay}")
                print(f"\t\tNum rows {p.nrow}")
                print(f"\t\tNum cols {p.ncol}")
                print(f"\t\tMax row res {np.asarray(p.delc).max()}")
                print(f"\t\tMax col res {np.asarray(p.delr).max()}")
                print(f"\t\tMin row res {np.asarray(p.delc).min()}")
                print(f"\t\tMin col res {np.asarray(p.delr).min()}")

    def plot_grid(self):
        fig = plt.figure(figsize=(10, 10))
        for i, d in enumerate(self.lst_subdomain_names):
            ax = fig.add_subplot(i + 1, 1, i + 1, aspect='equal')
            mgp = self.lst_subdomain_lgr[i].parent.modelgrid
            mgc = self.lst_subdomain_lgr[i].child.modelgrid
            mgc.plot(ax=ax, color='r')
            mgp.plot(ax=ax, color='b')
            fig.show()


class NestedDomainSimulation:
    def __init__(self,
                 sim,
                 parent_model_name,
                 lst_subdomain_names,
                 lst_subdomain_lgr,
                 ):
        try:
            self.sim = sim
            self.lst_subdomain_names = lst_subdomain_names
            self.lst_subdomain_lgr = lst_subdomain_lgr
            self.parent_model_name = parent_model_name

            self.streams_shp = None

            self._create_core_model_structure()

        except MFDataException as e:
            raise e

    def _setup_child_model(self, lgr, subdomain_name: str):
        logger.info(f"Creating core model structure for subdomain {subdomain_name}")
        lgrc = lgr.child
        gwf = flopy.mf6.ModflowGwf(self.sim, modelname=subdomain_name, save_flows=True)
        # Newton on with under relaxation
        gwf.name_file.newtonoptions = "UNDER_RELAXATION"
        dis = flopy.mf6.ModflowGwfdis(gwf, **lgrc.get_gridprops_dis6())
        oc = flopy.mf6.ModflowGwfoc(gwf,
                                    budget_filerecord=f"{subdomain_name}.cbb",
                                    head_filerecord=f"{subdomain_name}.hds",
                                    saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])

    def _create_core_model_structure(self):
        for lgr, name in zip(self.lst_subdomain_lgr, self.lst_subdomain_names):
            self._setup_child_model(lgr, name)

        try:
            ims_exists = isinstance(self.sim.ims, flopy.mf6.ModflowIms)
        except AttributeError:
            ims_exists = False

        if ims_exists:
            self.sim.ims.linear_acceleration = "bicgstab"  # self.sim.ims.build_mfdata("linear_acceleration", "bicgstab")
        else:
            ims_flow = flopy.mf6.ModflowIms(
                self.sim, linear_acceleration="BICGSTAB",
            )

        self.sim.register_ims_package(self.sim.ims, [self.parent_model_name] + self.lst_subdomain_names)

    def refine_grid_data(self,
                         streams_shp: str = None):
        self.streams_shp = streams_shp
        parent_model = self.sim.get_model(self.parent_model_name)
        for name, lgr in zip(self.lst_subdomain_names, self.lst_subdomain_lgr):
            for pck_name in parent_model.package_names:
                pck = parent_model.get_package(pck_name)
                if pck.package_type in ['ic', 'sto', 'npf', 'rcha', 'rch', 'maw', 'sfr']:
                    print(f"Package {pck_name} with {pck.package_type} found and will be regridded")
                    if pck.package_type == 'sfr' and streams_shp is None:
                        print("The parent model contains an SFR package but 'stream_shp' was not provided")
                        print("Skipping SFR regridding...")
                        continue
                    self.regrid_package(pck, parent_model, self.sim.get_model(name), lgr, name)

    @staticmethod
    def _regrid_data_layers(array3d, lgr) -> np.array:

        if array3d is None:
            return

        number_sublayers = lgr.ncppl  # if number of sublayers per parent layer is a constant
        lst_regridded_array = []
        for nlayp in range(lgr.nlayp):
            data_ = lgr.get_replicated_parent_array(array3d[nlayp, :, :])
            try:  # try if the number of sublayers per parent layers is a list
                number_sublayers = lgr.ncppl[nlayp]
            except:
                pass
            for _ in range(number_sublayers):
                lst_regridded_array.append(data_)

        cstrt = np.stack(lst_regridded_array)
        assert cstrt.shape == (lgr.nlay, lgr.nrow, lgr.ncol)
        return cstrt

    @staticmethod
    def _regrid_transient_layers(array3d, lgr) -> np.array:

        if array3d is None:
            return

        dct_regridded_array = {}
        for k, v in array3d.items():
            print(f"\tRegridding sp {k}")
            data_ = lgr.get_replicated_parent_array(v)
            dct_regridded_array[k] = data_

        return dct_regridded_array

    @singledispatchmethod
    def regrid_package(self, pmodel, cmodel, package, lgr, cname):
        raise NotImplementedError(package)

    @regrid_package.register
    def _(self,
          pkg: flopy.mf6.modflow.mfgwfic.ModflowGwfic,
          pmodel: flopy.mf6.ModflowGwf,
          cmodel: flopy.mf6.ModflowGwf,
          lgr: flopy.utils.lgrutil.Lgr,
          cname: str) -> flopy.mf6.modflow.mfgwfic.ModflowGwfic:

        strtp = pkg.strt.get_data()
        strtc = self._regrid_data_layers(strtp, lgr)

        cpkg = pkg.__class__(cmodel, strt=strtc)

        return cpkg

    @regrid_package.register
    def _(self,
          pkg: flopy.mf6.modflow.mfgwfsto.ModflowGwfsto,
          pmodel: flopy.mf6.modflow.ModflowGwf,
          cmodel: flopy.mf6.modflow.ModflowGwf,
          lgr: flopy.utils.lgrutil.Lgr,
          cname: str) -> flopy.mf6.modflow.mfgwfsto.ModflowGwfsto:

        data_arrays = [getattr(pkg, n).get_data() for n in
                       ["ss", "sy", "iconvert", ]]
        ssc, syc, iconvertc = [self._regrid_data_layers(data, lgr) for data
                               in data_arrays]

        if pkg.has_stress_period_data:
            raise NotImplementedError("regridding of data for tvs package not implemented")

        cpkg = pkg.__class__(cmodel,
                             storagecoefficient=pkg.storagecoefficient,
                             ss_confined_only=pkg.ss_confined_only,
                             steady_state=pkg.steady_state,
                             transient=pkg.transient,
                             ss=ssc,
                             sy=syc,
                             iconvert=iconvertc,
                             save_flows=True)

        return cpkg

    @regrid_package.register
    def _(self,
          pkg: flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf,
          pmodel: flopy.mf6.ModflowGwf,
          cmodel: flopy.mf6.ModflowGwf,
          lgr: flopy.utils.lgrutil.Lgr,
          cname: str) -> flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf:

        data_arrays = [getattr(pkg, n).get_data() for n in ["icelltype", "k", "k22", "k33",
                                                            "angle1", "angle2", "angle3", "wetdry"]]
        icelltypec, kc, k22c, k33c, angle1c, angle2c, angle3c, wetdryc = [self._regrid_data_layers(data, lgr)
                                                                          for data in data_arrays]

        cpkg = pkg.__class__(cmodel,
                             save_flows=True,
                             save_specific_discharge=pkg.save_specific_discharge,
                             alternative_cell_averaging=pkg.alternative_cell_averaging.data,
                             thickstrt=pkg.thickstrt.data,
                             cvoptions=pkg.cvoptions.data,
                             perched=pkg.perched.data,
                             k22overk=pkg.k22overk.data,
                             k33overk=pkg.k33overk.data,
                             dev_no_newton=pkg.dev_no_newton.data,
                             icelltype=icelltypec,
                             k=kc,
                             k22=k22c,
                             k33=k33c,
                             angle1=angle1c,
                             angle2=angle2c,
                             angle3=angle3c,
                             wetdry=wetdryc
                             )

        return cpkg

    @regrid_package.register
    def _(self,
          pkg: flopy.mf6.modflow.mfgwfrcha.ModflowGwfrcha,
          pmodel: flopy.mf6.ModflowGwf,
          cmodel: flopy.mf6.ModflowGwf,
          lgr: flopy.utils.lgrutil.Lgr,
          cname: str) -> flopy.mf6.modflow.mfgwfrcha.ModflowGwfrcha:

        print("About to regrid recharge transient grid...")
        rch_arrays = pkg.recharge.get_data()
        dct_rch = self._regrid_transient_layers(rch_arrays, lgr)

        cpkg = pkg.__class__(cmodel,
                             readasarrays=pkg.readasarrays,
                             recharge=dct_rch,
                             save_flows=True,
                             )

        return cpkg

    @regrid_package.register
    def _(self,
          pkg: flopy.mf6.modflow.mfgwfrch.ModflowGwfrch,
          pmodel: flopy.mf6.ModflowGwf,
          cmodel: flopy.mf6.ModflowGwf,
          lgr: flopy.utils.lgrutil.Lgr,
          cname: str) -> flopy.mf6.modflow.mfgwfrch.ModflowGwfrch:

        print("About to process recharge transient data...")
        rch_rec = pkg.stress_period_data
        lst_rch_rec = self._remap_stress_periods(rch_rec, lgr)

        cpkg = pkg.__class__(cmodel,
                             fixed_cell=pkg.fixed_cell,
                             auxiliary=pkg.auxiliary.array,
                             auxmultname=pkg.auxmultname.data,
                             stress_period_data=lst_rch_rec,
                             save_flows=True,
                             )
        print(cpkg.check())
        return cpkg

    @regrid_package.register
    def _regrid_package(self,
                        pkg: flopy.mf6.modflow.mfgwfmaw.ModflowGwfmaw,
                        pmodel: flopy.mf6.ModflowGwf,
                        cmodel: flopy.mf6.ModflowGwf,
                        lgr: flopy.utils.lgrutil.Lgr,
                        cname: str) -> flopy.mf6.modflow.mfgwfmaw.ModflowGwfmaw:

        fn_head_records = None
        fn_budget_records = None
        if pkg.head_filerecord.array is not None:
            fn_head_records = "cname_" + pkg.head_filerecord.array[0][0]
        if pkg.budget_filerecord.array is not None:
            fn_budget_records = "cname_" + pkg.budget_filerecord.array[0][0]

        print(f"Updating well connections in maw package {pkg.name}")
        # find subset of wells within the child grid
        lst_c_connections = []
        for w in pkg.connectiondata.array.T:
            l, r, c = w[2]
            if ((l >= lgr.nplbeg) &
                    (l <= lgr.nplend) &
                    (r >= lgr.nprbeg) &
                    (r <= lgr.nprend) &
                    (c >= lgr.npcbeg) &
                    (c <= lgr.npcend)):
                lst_c_connections.append(w)

        df_cconns = pd.DataFrame.from_records(lst_c_connections,
                                              columns=lst_c_connections[0].dtype.names)
        df_lst_rec = df_cconns.groupby(df_cconns['ifno']).apply(self._update_connection_records, lgr=lgr, include_groups=False)
        c_conns = df_lst_rec.to_records(index=False)

        # update ngwfnodes
        df_packagedata = pd.DataFrame.from_records(pkg.packagedata.array, columns=pkg.packagedata.dtype.names)
        df_c_packagedata = df_packagedata.loc[df_cconns['ifno'].unique()]
        df_c_packagedata['ngwfnodes'] = df_lst_rec.groupby(level=0).size()
        c_packagedata = df_c_packagedata.to_records(index=False)

        # update stress period data
        dct_c_sp = {}
        for i, sp in enumerate(pkg.perioddata.array):
            dct_c_sp[i] = pd.DataFrame.from_records(sp,
                                                    columns=sp.dtype.names
                                                    ).set_index('ifno').loc[df_c_packagedata.index].to_records(
                index=True)

        cpkg = pkg.__class__(cmodel,
                             save_flows=True,
                             print_input=pkg.print_input,
                             print_head=pkg.print_head,
                             boundnames=pkg.boundnames,
                             mover=pkg.mover,
                             head_filerecord=fn_head_records,
                             budget_filerecord=fn_budget_records,
                             no_well_storage=pkg.no_well_storage,
                             flow_correction=pkg.flow_correction,
                             flowing_wells=pkg.flowing_wells,
                             packagedata=c_packagedata,
                             connectiondata=c_conns,
                             perioddata=dct_c_sp
                             )
        return cpkg

    @regrid_package.register
    def _(self,
          pkg: flopy.mf6.modflow.mfgwfsfr.ModflowGwfsfr,
          pmodel: flopy.mf6.ModflowGwf,
          cmodel: flopy.mf6.ModflowGwf,
          lgr: flopy.utils.lgrutil.Lgr,
          cname: str) -> flopy.mf6.modflow.mfgwfsfr.ModflowGwfsfr:

        fn_stage_records = None
        fn_budget_records = None
        if pkg.stage_filerecord.array is not None:
            fn_stage_records = "cname_" + pkg.stage_filerecord.array[0][0]
        if pkg.budget_filerecord.array is not None:
            fn_budget_records = "cname_" + pkg.budget_filerecord.array[0][0]

        #reach_data = self.sfr_reach_data(pkg, lgr)

        df_conns = self._srf_develop_child_connections(pkg, lgr)

        pmodel_name = self.parent_model_name
        cmodel_name = cname,
        ppkg_name = pmodel.sfr.name
        cpkg_name = pkg.name

        # Retrieve mover package information before removing references to outside reaches and producing connection data
        mvr_period_data = self._generate_mover_period_data(pmodel_name, cmodel_name, ppkg_name, cpkg_name, df_conns)
        conn_data, ncons, ndvis = self._produce_package_connection_data(df_conns)

        package_data = self._produce_package_data(pkg, df_conns, ncons, ndvis)

        print("Instantiating new mover package to connect parent-child stream network")
        mvr_packages = [[pmodel_name, ppkg_name[0]], [cmodel_name[0], cpkg_name[0]]]
        mvr = self.new_mover_package(self.sim.gwfgwf,
                               maxmvr=max([len(per) for per in mvr_period_data.values()]),
                               packages=mvr_packages,
                               period_data=mvr_period_data,
                               name=f"{cname}_mvr.mvr"
                               )
        self.sim.gwfgwf.mvr_filerecord = f"{cname}_mvr.mvr"

        if not pmodel.sfr.mover:
            print(f"Updating parent model package {pmodel.sfr.name} with mover")
            pmodel.sfr.mover = True
            pmodel.sfr.write()

        cpkg = pkg.__class__(cmodel,
                             save_flows=pkg.save_flows,
                             stage_filerecord=fn_stage_records,
                             budget_filerecord=fn_budget_records,
                             mover=True,
                             maximum_picard_iterations=pkg.maximum_picard_iterations.data,
                             maximum_iterations=pkg.maximum_iterations.data,
                             maximum_depth_change=pkg.maximum_depth_change.data,
                             unit_conversion=pkg.unit_conversion.data,
                             length_conversion=pkg.length_conversion.data,
                             time_conversion=pkg.time_conversion.data,
                             connectiondata=conn_data,
                             packagedata=package_data,
                             perioddata={0: [[0, "INFLOW", 0.0]]}
                             )

        return cpkg

    @staticmethod
    def new_mover_package(model: flopy.mf6,
                          maxmvr: int,
                          packages: list,
                          period_data: list,
                          name: str) -> flopy.mf6.modflow.mfgwfmvr.ModflowGwfmvr:

        modelnames=False
        if isinstance(model, flopy.mf6.ModflowGwfgwf):
            modelnames = True

        cpkg = flopy.mf6.ModflowMvr(model,
                                    modelnames=modelnames,
                                    maxmvr=maxmvr,
                                    maxpackages=len(packages),
                                    packages=packages,
                                    perioddata=period_data,
                                    filename=f"{name}.mvr"
                                    )

        return cpkg


    @staticmethod
    def _produce_package_data(pkg, df_conns, ncons, ndivs):
        # Open SFR reach properties from parent model
        df_reach_data_p = pd.DataFrame.from_records(pkg.packagedata.array)

        # Merge with parent
        df_reach_data_c = df_conns.merge(df_reach_data_p, left_on='ifnop', right_on='ifno', how='left')

        # add layer number to the tuple of child cellid from parent cellids
        df_reach_data_c['cellids'] = [(df_reach_data_c.loc[i, 'cellid'][0],  # layer from parent cellid
                                       *df_reach_data_c.loc[i, 'cellids']) for i, _ in df_reach_data_c.iterrows()]

        # update number of connections, diversions, and fraction of incoming upstream flow
        df_reach_data_c['ncons'] = ncons
        df_reach_data_c['ndivs'] = ndivs

        rec_reach_data = df_reach_data_c[['cellids',
                         'lengths',
                         'rwid',
                         'rgrd',
                         'rtp',
                         'rbth',
                         'rhk',
                         'man',
                         'ncon',
                         'ustrf',
                         'ndv']].to_records(index=True)

        return rec_reach_data

    @staticmethod
    def _generate_mover_period_data(pmodel_name, cmodel_name, ppkg_name, cpkg_name, df_conns):

        absint = lambda x: int(np.abs(x))

        df_mover_conns = df_conns.filter(like='ic_').map(
            lambda x: x if abs(x) not in df_conns.index else np.nan
        ).dropna(
            axis=0,
            how='all'
        )

        parent_model_name = pmodel_name
        package_namep = ppkg_name[0]

        child_model_name = cmodel_name[0]
        package_namec = cpkg_name[0]

        mv_period_data = []
        for ifnoc, row in df_mover_conns.iterrows():
            for ifno in row:
                if math.isnan(ifno):
                    continue
                # Positive ifno indicates incoming water from parent
                if ifno >= 0:
                    mv_period_data.append(
                        [parent_model_name, package_namep, absint(ifno),
                         child_model_name, package_namec, absint(ifnoc),
                         "FACTOR", 1]
                    )
                else:
                    mv_period_data.append(
                        [child_model_name, package_namec, absint(ifnoc),
                         parent_model_name, package_namep, absint(ifno),
                         "FACTOR", 1]
                    )
        return {0: mv_period_data}

    @staticmethod
    def _produce_package_connection_data(df_conns):
        """
        Parses the connections dataframe, removes connections to reaches outside the child domain,
         removes NaNs and returns a list of tuples
        """
        # Apply to columns, drop NaNs from Series, and creates a list of tuples with network connections
        ls_conn_data = df_conns.filter(like='ic_').map(
            lambda x: x if abs(x) in df_conns.index else np.nan
        ).reset_index().T.apply(lambda x: tuple(x.dropna().astype(int))).to_list()

        nconns = [len(x)-1 for x in ls_conn_data]
        ndivs = [len(list(filter(lambda i: i<0, x)))-1 for x in ls_conn_data]

        return ls_conn_data, nconns, ndivs

    def _srf_develop_child_connections(self, pkg, lgr):

        # Open SFR connections from parent model
        df_conn = pd.DataFrame(pkg.connectiondata.array).set_index('ifno')

        # Open SFR reach properties from parent model
        df_reach_data_p = pd.DataFrame.from_records(pkg.packagedata.array)

        # get parent idomain
        idomainp = lgr.parent.idomain

        # retrieve child modelgrid and initialize a GridIntersect object
        mgrid = lgr.child.modelgrid
        ix = GridIntersect(mgrid, method='structured')

        df_shp_properties, ls_geoms = self._parse_shp_geom(self.streams_shp)

        # STAGE 1: develop connections between end child reaches that correspond to ends of parent reaches
        df_reach_conns = self._develop_segment_connections(df_shp_properties,
                                                           ls_geoms,
                                                           ix,
                                                           df_conn,
                                                           df_reach_data_p,
                                                           idomainp)

        # STAGE 2: Update reach ids and develop internal connections between segments of child reaches
        df_reach_conns = self._develop_internal_segment_connections(df_reach_conns)

        return df_reach_conns

    def _parse_shp_geom(self, fn_shp):
        """
        Open a vector dataset with the model river network
        Returns a dataframe with segment properties and a list of shapely geometries

        :param fn_shp: filename of vector dataset with river network
        :return: dataframe with reach properties, list of geometries (tuple)
        """

        with fiona.open(fn_shp) as shp_features:
            # Parse stream network properties as dataframe
            df_reach_data_shape_p = pd.DataFrame.from_records([feat.properties for feat in shp_features])
            ls_cgeom = [feat.geometry for feat in shp_features]

            return df_reach_data_shape_p, ls_cgeom

    @staticmethod
    def _develop_segment_connections(df_shp_properties, ls_geoms, ix, df_p_conn, df_p_props, idomainp):
        """
        Connects the end of child segments that correspond to the parent reaches and
        flags reach connections of high resolution reaches as NaNs

        :return: dataframe with child network connections based on parent connections
        """

        lst_cgeom = []
        for i, reachp in df_p_props.iterrows():
            ifno = reachp['ifno']
            idx = df_shp_properties[
                df_shp_properties['ReachID'] == ifno + 1].index.item()  # TODO: column names of properties
            cgeom = ix.intersect(shape(ls_geoms[idx]))
            cellidp = reachp['cellid']

            if (cgeom.size == 0) or (idomainp[cellidp] == 1):
                continue

            tmp_df = pd.DataFrame.from_records(cgeom)
            tmp_df['ifnop'] = ifno

            tmp_df = tmp_df.join(df_p_conn, on='ifnop')
            if tmp_df.shape[0] > 1:
                ls_col_idx = [tmp_df.T.index.get_loc(name) for name in tmp_df.iloc[0].filter(like='ic_').index]
                tmp_df.iloc[0, ls_col_idx] = tmp_df.iloc[0, ls_col_idx].map(lambda x: -REACH_FLAG if x < 0 else x)
                tmp_df.iloc[1:-1, ls_col_idx] = tmp_df.iloc[1:-1, ls_col_idx].map(lambda x: np.nan)
                tmp_df.iloc[-1, ls_col_idx] = tmp_df.iloc[-1, ls_col_idx].map(lambda x: REACH_FLAG if x > 0 else x)

            lst_cgeom.append(tmp_df)

        df_reach_data_c = pd.concat(lst_cgeom)
        df_reach_data_c = df_reach_data_c.reset_index(drop=True)
        df_reach_data_c.index.name = 'ifno'

        return df_reach_data_c

    def _develop_internal_segment_connections(self, df_reach_data_c):
        lst_cgeom = []
        for ifnop, tmp_df in df_reach_data_c.groupby('ifnop'):
            ls_col_idx = [tmp_df.T.index.get_loc(name) for name in tmp_df.iloc[0].filter(like='ic_').index]
            tmp_df.iloc[:, ls_col_idx] = tmp_df.iloc[:, ls_col_idx].map(self._map_connections_with_parent_segments,
                                                                        df=df_reach_data_c, na_action='ignore')
            if tmp_df.shape[0] > 1:
                tmp_df.iloc[0, ls_col_idx] = tmp_df.iloc[:1, ls_col_idx].map(
                    lambda x: -(tmp_df.index[0] + 1) if abs(x) == REACH_FLAG else x)
                tmp_df.iloc[1:-1, ls_col_idx[0]] = -(tmp_df.iloc[1:-1, ls_col_idx[0]].index + 1)
                tmp_df.iloc[1:-1, ls_col_idx[1]] = tmp_df.iloc[1:-1, ls_col_idx[1]].index - 1
                tmp_df.iloc[-1, ls_col_idx] = tmp_df.iloc[-1:, ls_col_idx].map(
                    lambda x: tmp_df.index[-1] - 1 if abs(x) == REACH_FLAG else x)

            lst_cgeom.append(tmp_df)

        return pd.concat(lst_cgeom)

    @staticmethod
    def _map_connections_with_parent_segments(x, df):
        df_subset = df[df['ifnop'] == abs(x)]
        if np.isnan(x) or df_subset.size == 0:
            return x
        if x < 0:
            return math.copysign(df_subset.index.min(), x)
        else:
            return math.copysign(df_subset.index.max(), x)

    def _remap_stress_periods(self, rch_rec, lgr):
        dct_recs = {}
        nlc, nrc, ncc = lgr.child.idomain.shape
        for sp, sp_data in rch_rec.data.items():
            lst_sp_rec = []
            for rec in sp_data:
                kp, ip, jp = rec[0]
                kc, ics, jcs = self.get_child_kij_index_connections(kp, ip, jp, lgr)
                for idx in range(len(ics)):
                    if any([(i<0 or i>=nrc or i>=ncc) for i in (kc[0], ics[idx], jcs[idx])]) or (lgr.child.idomain[kc[0], ics[idx], jcs[idx]] == 0):
                        continue
                    lst_sp_rec.append(
                        ((kc[0], ics[idx], jcs[idx]), *tuple(rec)[1:])
                    )
            dct_recs[sp] = lst_sp_rec
        return dct_recs

    @staticmethod
    def get_child_ij_indices(ip, jp, lgr):
        ic = (ip - lgr.nprbeg) * lgr.ncpp
        jc = (jp - lgr.npcbeg) * lgr.ncpp
        return ic, jc

    def get_child_kij_index_connections(self, kp, ip, jp, lgr):

        _, _, kcs = self.get_child_layer_connections(0, kp, lgr)
        ic = (ip - lgr.nprbeg) * lgr.ncpp
        jc = (jp - lgr.npcbeg) * lgr.ncpp
        return kcs, np.arange(ic, ic+lgr.ncpp), np.arange(jc, jc+lgr.ncpp)

    @staticmethod
    def get_child_layer_connections(icounter, kp, lgr):
        n_sublayers = lgr.ncppl[kp]
        iconns = np.arange(icounter, icounter + n_sublayers)
        icounter += n_sublayers
        lstart = lgr.ncppl.cumsum()[kp] - n_sublayers
        kcs = np.arange(lstart, lgr.ncppl.cumsum()[kp])
        return icounter, iconns, kcs

    def _update_connection_records(self, recs, lgr):
        next_iconn = 0
        lst_recs = []
        for i, rec in recs.iterrows():

            kp, ip, jp = rec['cellid']
            ic, jc = self.get_child_ij_indices(ip, jp, lgr)
            next_iconn, iconns, clyrs = self.get_child_layer_connections(next_iconn, kp, lgr)
            for iconn, kc in zip(iconns, clyrs):
                r = np.copy(rec)
                # copied array is not longer a rec array.
                # Assign values by position index
                r[1] = iconn
                r[2] = (kc, ic, jc)
                lst_recs.append(r)

        return pd.DataFrame(lst_recs)

    def write_simulation(self, sim_ws: str = None):
        if sim_ws is not None:
            self.sim.set_sim_path(sim_ws)

        self.sim.write_simulation()
