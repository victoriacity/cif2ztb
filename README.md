# cif2ztb
Tabulated zeolite potential

readcif: convert CIF file to cell info and coordinates
* coordinates: uint8[] atom_id, float[] x_coord, float[] y_coord, float[] z_coord
* cell info: float[3] cell_param, float[3] cell_angle

calc_energy: use CUDA to calculate tabulated energies

write_h5: output tabulated potential to HDF5 file

write_ztb: output tabulated potential to ZTB file
