import unittest
import matplotlib.pyplot as plt
from src.eki_mf6_utils import NestedDomain, NestedDomainSimulation
import pickle
import copy

from flopy.mf6.modflow.mfgwfgwf import ModflowGwfgwf


class TestNestedDomain(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.nested_grid = NestedDomain.from_parent_model(sim_ws="./data/gwf")
        with open("./data/gwf/lgr.pckl", "rb") as file:
            cls.lgr = pickle.load(file)

    def test_constructor_filenotfound(self):
        with self.assertRaises(FileNotFoundError):
            NestedDomain.from_parent_model()

    def test_constructor(self):
        sim = NestedDomain.from_parent_model(sim_ws="./data/gwf")
        assert isinstance(sim, NestedDomain)

    def test_define_subdomain(self):
        # properties of model domain
        xoff = 6135367
        yoff = 2047406
        crs = "EPSG:2227"
        angrot = 0

        kstart = 0
        kstop = 10
        istart = 60
        istop = 70
        jstart = 50
        jstop = 60
        self.nested_grid.define_subdomain(name="subdomain",
                                          kstart=kstart,
                                          kstop=kstop,
                                          istart=istart,
                                          istop=istop,
                                          jstart=jstart,
                                          jstop=jstop,
                                          xoff=xoff,
                                          yoff=yoff,
                                          angrot=angrot,
                                          crs=crs,
                                          num_cells_per_parent_cell=2,
                                          num_layers_per_parent_layer=[2, 2, 2] + [1]*7)
        assert (self.nested_grid.gwf.dis.idomain.get_data()[kstart:kstop+1, istart:istop+1, jstart:jstop+1] == 0).all()
        assert self.nested_grid.lst_subdomain_lgr[0] is not None

    def test_define_subdomain_shp(self):
        shp ='data/gwf/subdomain.shp'
        self.nested_grid.define_subdomain(name="subdomain",
                                          nested_domain_shp=shp,
                                          xoff=6135367,
                                          yoff=2047406,
                                          num_cells_per_parent_cell=2,
                                          num_layers_per_parent_layer=[2, 2, 2] + [1]*7)
        plt.imshow(self.nested_grid.gwf.dis.idomain.get_data()[1,:,:])
        plt.show()
        assert self.nested_grid.lst_subdomain_lgr[0] is not None

    def test_plot_grid(self):
        nd = NestedDomain()
        nd.lst_subdomain_names.append("some_subdomain")
        nd.lst_subdomain_lgr.append(self.lgr)
        nd.plot_grid()

    def test_create_exchange_data(self):
        nd = NestedDomain()
        nd.gwf = self.nested_grid.gwf
        nd.sim = copy.deepcopy(self.nested_grid.sim)
        ret = nd._create_exchange_data(self.lgr, "test_subdomain")
        assert isinstance(ret, ModflowGwfgwf)

    def test_nested_grid_happy_path(self):
        kstart = 0
        kstop = 10
        istart = 60
        istop = 70
        jstart = 50
        jstop = 60
        self.nested_grid.define_subdomain(name="subdomain",
                                  kstart=kstart,
                                  kstop=kstop,
                                  istart=istart,
                                  istop=istop,
                                  jstart=jstart,
                                  jstop=jstop,
                                  num_cells_per_parent_cell=2,
                                  num_layers_per_parent_layer=[2, 2, 2] + [1]*7)
        ng_sim = self.nested_grid.get_flow_simulation()
        assert isinstance(ng_sim, NestedDomainSimulation)
        assert ng_sim.sim is not None
        assert len(ng_sim.lst_subdomain_names) > 0
        assert len(ng_sim.lst_subdomain_lgr) > 0
        assert ("gwfgwf" in [pkg.package_type for pkg in ng_sim.sim.sim_package_list])


class TestNestedDomainSimulation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        with open("./data/gwf/ng_simul.pckl", "rb") as file:
            cls.nd_sim = pickle.load(file)

    def test_refine_grid_data(self):
        self.nd_sim.refine_grid_data()

    def test_refine_grid_data_with_sfr(self):

        with open('data/gwf/nd_simul_irregu_domain_sfr.pckl', 'rb') as f:
            nd_sim = pickle.load(f)

        # # load the sfr package into the unpickled model
        # nested_domain = NestedDomain.from_parent_model(sim_ws="./data/gwf")
        # # self.nd_sim.sim.get_model().load_package(ftype='sfr', fname='Zone7_gwm_2024.sfr', pname='sfr',
        # #                                         strict=True, ref_path='.')
        # # self.nd_sim.streams_shp = './data/gwf/SFR_reaches_final_v2.shp'
        #
        #
        # shp = 'data/gwf/subdomain.shp'
        # nested_domain.define_subdomain(name="subdomain",
        #                                nested_domain_shp=shp,
        #                                xoff=6135367,
        #                                yoff=2047406,
        #                                num_cells_per_parent_cell=2,
        #                                num_layers_per_parent_layer=[2, 2, 2] + [1] * 7)
        # nd_sim = nested_domain.get_flow_simulation()

        nd_sim.refine_grid_data(streams_shp='./data/gwf/SFR_reaches_final_v2.shp')

    def test_write_simulation(self):
        self.nd_sim.write_simulation(sim_ws='./data/test_gwf')


if __name__ == '__main__':
    unittest.main()

