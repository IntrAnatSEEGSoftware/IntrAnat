from soma import aims
import pdb

mesh_filename = '/home/b67-belledone/Desktop/epilepsie-manik/MNI_Brainvisa/t1mri/T1pre_1900-1-3/default_analysis/segmentation/mesh/inflatedGre_2016_MNI1_Lwhite.gii'
regions_texture_filename = '/home/b67-belledone/Desktop/epilepsie-manik/MNI_Brainvisa/t1mri/T1pre_1900-1-3/default_analysis/segmentation/mesh/surface_analysis/Gre_2016_MNI1_Lwhite_parcels_marsAtlas.gii'
boundaries_filename = '/tmp/inflated_marsatlas_boundariesL.mesh'

mesh  = aims.read(mesh_filename)
regions_texture = aims.read(regions_texture_filename)
boundaries = aims.SurfaceManip.meshTextureBoundary(mesh, regions_texture, -1)
boundaries.header()['material'] =  { 'diffuse': [1., 0., 0., 1.], 'line_width': 6}
aims.write(boundaries, boundaries_filename)