
EKI_MF6_UTILS
=============

Collection of utilities to manipulate and facilitate the development of MODFLOW 6 models. 

Currently, the main functionality implemented in the library are two classes  designed to enable nested domain modeling of groundwater system using MODFLOW 6.

* NestedDomain
* NestedDomainSimulation 

These two libraries allow the reation of high-resolution subdomains within a coarser parent model grid, enabling detailed local modeling while maintaining computational efficiency.

As of now, the libraries uses an internal version of `flopy` to create nested subdomains with arbitrary shapes. The library maps grid connections and has functionality
to regrid the following MODFLOW package types: IC, STO, NPF, RCH, MAW, SFR, and CHD


USAGE
=====

Installation
------------

Create a new conda environment:

```conda create -n eki_mf6_utils python=3.11```

Activate the environment 

```conda activate eki_mf6_utils```

Install the library

```pip install git+https://github.com/mmaneta/eki_mf6_utils.git@develop ```


Example 
--------

In a python notebook import the required libraries

```angular2html
#import flopy.mf6.mfsimbase
#%matplotlib nbagg
#import numpy as np
import eki_lgr_utils.nested_domain as ekind
#import mpld3
```

Create the base of the nested domain from the parent MF6 model:   

```
model_dir=r"Z:\Zone7WaterAgency_C20037\WorkEfforts\Model\LocalGridRefinement\test_lgr_model\gwf_PackageTest11"
nested_simulation = ekind.NestedDomain.from_parent_model(sim_ws=model_dir)
```

There are two ways to designate the area where the grid will be refined: 

1. For rectangular domains, information on the area to be refined can be bound
by the beginning and end row (i), columns (j) and layer (k). 

```angular2html
kstart = 0
kstop = 10
istart = 50
istop = 70
jstart = 35
jstop = 55
nested_simulation.define_subdomain(name="subdomain",
                                  kstart=kstart,
                                  kstop=kstop,
                                  istart=istart,
                                  istop=istop,
                                  jstart=jstart,
                                  jstop=jstop,
                                  num_cells_per_parent_cell=3,
                                  num_layers_per_parent_layer=1#[2, 2, 2] + [1]*7
                            )
```

2. Alternatively, for regions that are not rectangular, the domain can be defined with 
a vector polygon file (shapefile, geojson, among other compatible formats). The polygon must have a
a subdomain ID integer field (`feature_field`). The cooordinates of the lower left corner of the 
parent grid as well as the grid rotation need to be provided. The parent grid needs to be in the same projection (crs) 
as `nested_domain_shp`.


```angular2html
# properties of model domain
xoff = 6135367
yoff = 2033406
angrot = 0

fn_shapefile = r"Z:\Zone7WaterAgency_C20037\WorkEfforts\Model\LocalGridRefinement\test_lgr_model\shapefiles\pfos_area_polygon.shp"
nested_simulation.define_subdomain(name="pfos_subdomain",
                            xoff=xoff, 
                            yoff=yoff,                             
                            angrot=angrot,
                            nested_domain_shp=fn_shapefile,
                            feature_field='LGR_PFOS',
                            num_cells_per_parent_cell=3,
                            num_layers_per_parent_layer=1
                            )
        
```

Once the subdomain is defined, generate the basic flow model for the parent and nested domains along with flow exchange
information:


```angular2html
nested_simulation_flow_model = nested_simulation.get_flow_simulation()
```

```angular2html
fn_stream_shp=r"Z:\Zone7WaterAgency_C20037\WorkEfforts\Model\LocalGridRefinement\test_lgr_model\shapefiles\SFR_reaches_final_v2.shp"
nested_simulation_flow_model.refine_grid_data(streams_shp=fn_stream_shp)
```


A methods to generate a transport model will be implemented in the future

Once the model is generated you can write the model files to the hard drive:

```angular2html
nested_simulation_flow_model.write_simulation('./Zone7_PFOS_gwf')
```