from functools import singledispatchmethod
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import flopy
from flopy.mf6.mfbase import MFDataException
from flopy.utils.lgrutil import Lgr

from logging import getLogger

logger = getLogger(__name__)


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
                         istart: int,
                         istop: int,
                         jstart: int,
                         jstop: int,
                         kstart: int,
                         kstop: int,
                         num_cells_per_parent_cell: int = 3,
                         num_layers_per_parent_layer: list = 1):

        self.lst_subdomain_names.append(name)
        idomain = self.gwf.dis.idomain.get_data()

        ## deactivate the cells in the parent domain where the child grid will be placed
        idomain[kstart:kstop + 1, istart:istop + 1, jstart:jstop + 1] = 0
        self.gwf.dis.idomain.set_data(idomain)
        idomain[kstart:kstop + 1, istart:istop + 1, jstart:jstop + 1] = 2

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
                                      exchangedata=exchangedata
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
        for i,d in enumerate(self.lst_subdomain_names):

            print(f"DOMAIN {d}:")
            for name, p in zip(["Parent grid", "Nested grid"], [self.lst_subdomain_lgr[i].parent, self.lst_subdomain_lgr[i].child]):

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
        for i,d in enumerate(self.lst_subdomain_names):
            ax = fig.add_subplot(i+1, 1, i+1, aspect='equal')
            mgp = self.lst_subdomain_lgr[i].parent.modelgrid
            mgc = self.lst_subdomain_lgr[i].child.modelgrid
            mgc.plot(ax=ax, color='r')
            mgp.plot(ax=ax, color='b')
            fig.show()


class NestedDomainSimulation:
    def __init__(self, sim, parent_model_name, lst_subdomain_names, lst_subdomain_lgr):
        try:
            self.sim = sim
            self.lst_subdomain_names = lst_subdomain_names
            self.lst_subdomain_lgr = lst_subdomain_lgr
            self.parent_model_name = parent_model_name

            self._create_core_model_structure()

        except MFDataException as e:
            raise e

    def _setup_child_model(self, lgr, subdomain_name: str):
        logger.info(f"Creating core model structure for subdomain {subdomain_name}")
        lgrc = lgr.child
        gwf = flopy.mf6.ModflowGwf(self.sim, modelname=subdomain_name, save_flows=True)
        #newton on with under relaxation
        gwf.name_file.newtonoptions="UNDER_RELAXATION"
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
            self.sim.ims.linear_acceleration = "bicgstab"  #self.sim.ims.build_mfdata("linear_acceleration", "bicgstab")
        else:
            ims_flow = flopy.mf6.ModflowIms(
                self.sim, linear_acceleration="BICGSTAB",
            )

        self.sim.register_ims_package(self.sim.ims, [self.parent_model_name] + self.lst_subdomain_names)

    def refine_grid_data(self):
        parent_model = self.sim.get_model(self.parent_model_name)
        for name, lgr in zip(self.lst_subdomain_names, self.lst_subdomain_lgr):
            for pck_name in parent_model.package_names:
                pck = parent_model.get_package(pck_name)
                if pck.package_type in ['ic', 'sto', 'npf', 'rcha', 'maw']:
                    logger.info(f"Package {pck_name} with {pck.package_type} found and will be regridded")
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
        for k,v in array3d.items():
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
    def _regrid_package(self,
                        pkg: flopy.mf6.modflow.mfgwfmaw.ModflowGwfmaw,
                        pmodel: flopy.mf6.ModflowGwf,
                        cmodel: flopy.mf6.ModflowGwf,
                        lgr: flopy.utils.lgrutil.Lgr,
                        cname: str) -> flopy.mf6.modflow.mfgwfmaw.ModflowGwfmaw:

        fn_head_records = None
        fn_budget_records = None
        if pkg.head_filerecord.array is not None:
            fn_head_records = "cname_" + pkg.head_filerecord
        if pkg.budget_filerecord.array is not None:
            fn_budget_records = "cname_" + pkg.budget_filerecord

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
        df_lst_rec = df_cconns.groupby(df_cconns['ifno']).apply(self._update_connection_records, lgr=lgr)
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
                                          ).set_index('ifno').loc[df_c_packagedata.index].to_records(index=True)


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

    @staticmethod
    def get_child_ij_indices(ip, jp, lgr):
        ic = (ip - lgr.nprbeg) * lgr.ncpp
        jc = (jp - lgr.npcbeg) * lgr.ncpp
        return ic, jc

    @staticmethod
    def get_child_layer_connections(icounter, kp, lgr):
        n_sublayers = lgr.ncppl[kp]
        iconns = np.arange(icounter, icounter+n_sublayers)
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







