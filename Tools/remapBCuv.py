import os
import numpy as np
from netCDF4 import Dataset, date2num, num2date

import pyroms
import pyroms_toolbox

class nctime(object):
    pass

def remapBCuv(src_fileuv, src_vnameuv, wts_file, src_grd, dst_grd, dst_fileuv, dxy=20, cdepth=0, kk=2, verbose=True):

    src_ncuv = Dataset(src_fileuv)
    src_varu = src_ncuv.variables[src_vnameuv[0]][0]
    src_varv = src_ncuv.variables[src_vnameuv[1]][0]
    # read time and change units to 'days' for ROMS
    dtime = num2date(src_ncuv.variables['time'][0],units=src_ncuv['time'].units,calendar=src_ncuv['time'].calendar)
    time  = date2num(dtime,units='days since 2000-01-01 00:00:00',calendar='gregorian')
    
    # get time
    nctime.long_name = 'time'
    nctime.units = 'days' # for ROMS
    #nctime.calendar = 'gregorian'

    # get dimensions
    Mp, Lp = dst_grd.hgrid.mask_rho.shape

    # create destination file
    for f in dst_fileuv:
        print('==>Info: Creating destination file', f)
        if os.path.exists(f) is True: os.remove(f)
        pyroms_toolbox.nc_create_roms_file(f, dst_grd, nctime)

    # open destination file
    ncu = Dataset(dst_fileuv[0], 'a', format='NETCDF3_64BIT')
    ncv = Dataset(dst_fileuv[1], 'a', format='NETCDF3_64BIT')

    #get missing value
    spval = src_ncuv[src_vnameuv[0]]._FillValue
    src_varu = src_ncuv[src_vnameuv[0]][0,:,1:-1,1:-1] # remove time and 4 corners
    src_varv = src_ncuv[src_vnameuv[1]][0,:,1:-1,1:-1] # remove time and 4 corners

    # get weights file
    #wts_file = 'remap_weights_%s_to_%s_bilinear_t_to_rho.nc'%(src_name,dst_name)

    # build intermediate zgrid
    zlevel = -src_grd.z_t[::-1,0,0]
    nzlevel = len(zlevel)
    dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
    dst_grdz = pyroms.grid.ROMS_Grid(dst_fileuv[0]+'_Z', dst_grd.hgrid, dst_zcoord)
    
    # create variable in destination file

    if verbose: print('==>Info: Creating variable u_north')
    ncu.createVariable('u_north', 'f8', ('ocean_time', 's_rho', 'xi_u'), fill_value=spval)
    ncu.variables['u_north'].long_name = '3D u-momentum north boundary condition'
    ncu.variables['u_north'].units = 'meter second-1'
    ncu.variables['u_north'].field = 'u_north, scalar, series'
    if verbose: print('==>Info: Creating variable u_south')
    ncu.createVariable('u_south', 'f8', ('ocean_time', 's_rho', 'xi_u'), fill_value=spval)
    ncu.variables['u_south'].long_name = '3D u-momentum south boundary condition'
    ncu.variables['u_south'].units = 'meter second-1'
    ncu.variables['u_south'].field = 'u_south, scalar, series'
    if verbose: print('==>Info: Creating variable u_east')
    ncu.createVariable('u_east', 'f8', ('ocean_time', 's_rho', 'eta_u'), fill_value=spval)
    ncu.variables['u_east'].long_name = '3D u-momentum east boundary condition'
    ncu.variables['u_east'].units = 'meter second-1'
    ncu.variables['u_east'].field = 'u_east, scalar, series'
    if verbose: print('==>Info: Creating variable u_west')
    ncu.createVariable('u_west', 'f8', ('ocean_time', 's_rho', 'eta_u'), fill_value=spval)
    ncu.variables['u_west'].long_name = '3D u-momentum west boundary condition'
    ncu.variables['u_west'].units = 'meter second-1'
    ncu.variables['u_west'].field = 'u_east, scalar, series'
    
    if verbose: print('==>Info: Creating variable ubar_north')
    ncu.createVariable('ubar_north', 'f8', ('ocean_time', 'xi_u'), fill_value=spval)
    ncu.variables['ubar_north'].long_name = '2D u-momentum north boundary condition'
    ncu.variables['ubar_north'].units = 'meter second-1'
    ncu.variables['ubar_north'].field = 'ubar_north, scalar, series'
    if verbose: print('==>Info: Creating variable ubar_south')
    ncu.createVariable('ubar_south', 'f8', ('ocean_time', 'xi_u'), fill_value=spval)
    ncu.variables['ubar_south'].long_name = '2D u-momentum south boundary condition'
    ncu.variables['ubar_south'].units = 'meter second-1'
    ncu.variables['ubar_south'].field = 'ubar_south, scalar, series'
    if verbose: print('==>Info: Creating variable ubar_east')
    ncu.createVariable('ubar_east', 'f8', ('ocean_time', 'eta_u'), fill_value=spval)
    ncu.variables['ubar_east'].long_name = '2D u-momentum east boundary condition'
    ncu.variables['ubar_east'].units = 'meter second-1'
    ncu.variables['ubar_east'].field = 'ubar_east, scalar, series'
    if verbose: print('==>Info: Creating variable ubar_west')
    ncu.createVariable('ubar_west', 'f8', ('ocean_time', 'eta_u'), fill_value=spval)
    ncu.variables['ubar_west'].long_name = '2D u-momentum west boundary condition'
    ncu.variables['ubar_west'].units = 'meter second-1'
    ncu.variables['ubar_west'].field = 'ubar_east, scalar, series'
    
    if verbose: print('==>Info: Creating variable v_north')
    ncv.createVariable('v_north', 'f8', ('ocean_time', 's_rho', 'xi_v'), fill_value=spval)
    ncv.variables['v_north'].long_name = '3D v-momentum north boundary condition'
    ncv.variables['v_north'].units = 'meter second-1'
    ncv.variables['v_north'].field = 'v_north, scalar, series'
    if verbose: print('==>Info: Creating variable v_south')
    ncv.createVariable('v_south', 'f8', ('ocean_time', 's_rho', 'xi_v'), fill_value=spval)
    ncv.variables['v_south'].long_name = '3D v-momentum south boundary condition'
    ncv.variables['v_south'].units = 'meter second-1'
    ncv.variables['v_south'].field = 'v_south, scalar, series'
    if verbose: print('==>Info: Creating variable v_east')
    ncv.createVariable('v_east', 'f8', ('ocean_time', 's_rho', 'eta_v'), fill_value=spval)
    ncv.variables['v_east'].long_name = '3D v-momentum east boundary condition'
    ncv.variables['v_east'].units = 'meter second-1'
    ncv.variables['v_east'].field = 'v_east, scalar, series'
    if verbose: print('==>Info: Creating variable v_west')
    ncv.createVariable('v_west', 'f8', ('ocean_time', 's_rho', 'eta_v'), fill_value=spval)
    ncv.variables['v_west'].long_name = '3D v-momentum west boundary condition'
    ncv.variables['v_west'].units = 'meter second-1'
    ncv.variables['v_west'].field = 'v_east, scalar, series'
    
    if verbose: print('==>Info: Creating variable vbar_north')
    ncv.createVariable('vbar_north', 'f8', ('ocean_time', 'xi_v'), fill_value=spval)
    ncv.variables['vbar_north'].long_name = '2D v-momentum north boundary condition'
    ncv.variables['vbar_north'].units = 'meter second-1'
    ncv.variables['vbar_north'].field = 'vbar_north, scalar, series'
    if verbose: print('==>Info: Creating variable vbar_south')
    ncv.createVariable('vbar_south', 'f8', ('ocean_time', 'xi_v'), fill_value=spval)
    ncv.variables['vbar_south'].long_name = '2D v-momentum south boundary condition'
    ncv.variables['vbar_south'].units = 'meter second-1'
    ncv.variables['vbar_south'].field = 'vbar_south, scalar, series'
    if verbose: print('==>Info: Creating variable vbar_east')
    ncv.createVariable('vbar_east', 'f8', ('ocean_time', 'eta_v'), fill_value=spval)
    ncv.variables['vbar_east'].long_name = '2D v-momentum east boundary condition'
    ncv.variables['vbar_east'].units = 'meter second-1'
    ncv.variables['vbar_east'].field = 'vbar_east, scalar, series'
    if verbose: print('Creating variable vbar_west')
    ncv.createVariable('vbar_west', 'f8', ('ocean_time', 'eta_v'), fill_value=spval)
    ncv.variables['vbar_west'].long_name = '2D v-momentum west boundary condition'
    ncv.variables['vbar_west'].units = 'meter second-1'
    ncv.variables['vbar_west'].field = 'vbar_east, scalar, series'
    
    
    # remaping
    if verbose: print('==>Info: remapping and rotating u and v for B.C.', ', time =', time)


    # flood the grid
    if verbose: print('==>Info: flood the grid')
    src_uz = pyroms_toolbox.Grid_HYCOM.flood_fast(src_varu, src_grd, pos='t', \
                spval=spval, dxy=dxy, cdepth=cdepth, kk=kk, verbose=verbose)
    src_vz = pyroms_toolbox.Grid_HYCOM.flood_fast(src_varv, src_grd, pos='t', \
                spval=spval, dxy=dxy, cdepth=cdepth, kk=kk, verbose=verbose)

    # horizontal interpolation using scrip weights
    if verbose: print('==>Info: horizontal interpolation using scrip weights')
    dst_uz = pyroms.remapping.remap(src_uz, wts_file, spval=spval)
    dst_vz = pyroms.remapping.remap(src_vz, wts_file, spval=spval)

    # vertical interpolation from standard z level to sigma
    if verbose: print('==>Info: vertical interpolation from standard z level to sigma')
    dst_u_north = pyroms.remapping.z2roms(dst_uz[::-1, Mp-2:Mp, 0:Lp], \
                      dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                      flood=False, irange=(0,Lp), jrange=(Mp-2,Mp))
    dst_u_south = pyroms.remapping.z2roms(dst_uz[::-1, 0:2, 0:Lp], \
                      dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                      flood=False, irange=(0,Lp), jrange=(0,2))
    dst_u_east = pyroms.remapping.z2roms(dst_uz[::-1, 0:Mp, Lp-2:Lp], \
                      dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                      flood=False, irange=(Lp-2,Lp), jrange=(0,Mp))
    dst_u_west = pyroms.remapping.z2roms(dst_uz[::-1, 0:Mp, 0:2], \
                      dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                      flood=False, irange=(0,2), jrange=(0,Mp))
    
    dst_v_north = pyroms.remapping.z2roms(dst_vz[::-1, Mp-2:Mp, 0:Lp], \
                      dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                      flood=False, irange=(0,Lp), jrange=(Mp-2,Mp))
    dst_v_south = pyroms.remapping.z2roms(dst_vz[::-1, 0:2, 0:Lp], \
                      dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                      flood=False, irange=(0,Lp), jrange=(0,2))
    dst_v_east = pyroms.remapping.z2roms(dst_vz[::-1, 0:Mp, Lp-2:Lp], \
                      dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                      flood=False, irange=(Lp-2,Lp), jrange=(0,Mp))
    dst_v_west = pyroms.remapping.z2roms(dst_vz[::-1, 0:Mp, 0:2], \
                      dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                      flood=False, irange=(0,2), jrange=(0,Mp))

    # rotate u,v fields
    src_angle = pyroms.remapping.remap(src_grd.angle, wts_file, spval=spval)

    dst_angle = dst_grd.hgrid.angle_rho
    angle = dst_angle - src_angle
    angle = np.tile(angle, (dst_grd.vgrid.N, 1, 1))

    U_north = dst_u_north + dst_v_north*1j
    eitheta_north = np.exp(-1j*angle[:,Mp-2:Mp, 0:Lp])
    U_north = U_north * eitheta_north
    dst_u_north = np.real(U_north)
    dst_v_north = np.imag(U_north)

    U_south = dst_u_south + dst_v_south*1j
    eitheta_south = np.exp(-1j*angle[:,0:2, 0:Lp])
    U_south = U_south * eitheta_south
    dst_u_south = np.real(U_south)
    dst_v_south = np.imag(U_south)

    U_east = dst_u_east + dst_v_east*1j
    eitheta_east = np.exp(-1j*angle[:,0:Mp, Lp-2:Lp])
    U_east = U_east * eitheta_east
    dst_u_east = np.real(U_east)
    dst_v_east = np.imag(U_east)
    
    U_west = dst_u_west + dst_v_west*1j
    eitheta_west = np.exp(-1j*angle[:,0:Mp, 0:2])
    U_west = U_west * eitheta_west
    dst_u_west = np.real(U_west)
    dst_v_west = np.imag(U_west)
    
    # move back to u,v points
    dst_u_north = 0.5 * np.squeeze(dst_u_north[:,-1,:-1] + dst_u_north[:,-1,1:])
    dst_v_north = 0.5 * np.squeeze(dst_v_north[:,:-1,:] + dst_v_north[:,1:,:])
    dst_u_south = 0.5 * np.squeeze(dst_u_south[:,0,:-1] + dst_u_south[:,0,1:])
    dst_v_south = 0.5 * np.squeeze(dst_v_south[:,:-1,:] + dst_v_south[:,1:,:])
    dst_u_east = 0.5 * np.squeeze(dst_u_east[:,:,:-1] + dst_u_east[:,:,1:])
    dst_v_east = 0.5 * np.squeeze(dst_v_east[:,:-1,-1] + dst_v_east[:,1:,-1])
    dst_u_west = 0.5 * np.squeeze(dst_u_west[:,:,:-1] + dst_u_west[:,:,1:])
    dst_v_west = 0.5 * np.squeeze(dst_v_west[:,:-1,0] + dst_v_west[:,1:,0])

    # Fix north pole
#    print dst_u_north.shape
#    print dst_v_north.shape
#    dst_u_north[:,707] = 2/3.0*dst_u_north[:,706] + 1/3.0*dst_u_north[:,709]
#    dst_u_north[:,708] = 1/3.0*dst_u_north[:,706] + 2/3.0*dst_u_north[:,709]
#    dst_v_north[:,708] = 2/3.0*dst_v_north[:,707] + 1/3.0*dst_v_north[:,710]
#    dst_v_north[:,709] = 1/3.0*dst_v_north[:,707] + 2/3.0*dst_v_north[:,710]

    # spval
    idxu_north = np.where(dst_grd.hgrid.mask_u[-1,:] == 0)
    idxv_north = np.where(dst_grd.hgrid.mask_v[-1,:] == 0)
    idxu_south = np.where(dst_grd.hgrid.mask_u[0,:] == 0)
    idxv_south = np.where(dst_grd.hgrid.mask_v[0,:] == 0)
    idxu_east = np.where(dst_grd.hgrid.mask_u[:,-1] == 0)
    idxv_east = np.where(dst_grd.hgrid.mask_v[:,-1] == 0)
    idxu_west = np.where(dst_grd.hgrid.mask_u[:,0] == 0)
    idxv_west = np.where(dst_grd.hgrid.mask_v[:,0] == 0)
    for n in range(dst_grd.vgrid.N):
        dst_u_north[n, idxu_north[0]] = spval
        dst_v_north[n, idxv_north[0]] = spval
        dst_u_south[n, idxu_south[0]] = spval
        dst_v_south[n, idxv_south[0]] = spval
        dst_u_east[n, idxu_east[0]] = spval
        dst_v_east[n, idxv_east[0]] = spval
        dst_u_west[n, idxu_west[0]] = spval
        dst_v_west[n, idxv_west[0]] = spval

    # compute depth average velocity ubar and vbar
    # get z at the right position
    z_u_north = 0.5 * (dst_grd.vgrid.z_w[0,:,-1,:-1] + dst_grd.vgrid.z_w[0,:,-1,1:])
    z_v_north = 0.5 * (dst_grd.vgrid.z_w[0,:,-1,:] + dst_grd.vgrid.z_w[0,:,-2,:])
    z_u_south = 0.5 * (dst_grd.vgrid.z_w[0,:,0,:-1] + dst_grd.vgrid.z_w[0,:,0,1:])
    z_v_south = 0.5 * (dst_grd.vgrid.z_w[0,:,0,:] + dst_grd.vgrid.z_w[0,:,1,:])
    z_u_east = 0.5 * (dst_grd.vgrid.z_w[0,:,:,-1] + dst_grd.vgrid.z_w[0,:,:,-2])
    z_v_east = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,-1] + dst_grd.vgrid.z_w[0,:,1:,-1])
    z_u_west = 0.5 * (dst_grd.vgrid.z_w[0,:,:,0] + dst_grd.vgrid.z_w[0,:,:,1])
    z_v_west = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,0] + dst_grd.vgrid.z_w[0,:,1:,0])

    dst_ubar_north = np.zeros(dst_u_north.shape[1])
    dst_ubar_south = np.zeros(dst_u_south.shape[1])
    dst_ubar_east = np.zeros(dst_u_east.shape[1])
    dst_ubar_west = np.zeros(dst_u_west.shape[1])
    dst_vbar_north = np.zeros(dst_v_north.shape[1])
    dst_vbar_south = np.zeros(dst_v_south.shape[1])
    dst_vbar_east = np.zeros(dst_v_east.shape[1])
    dst_vbar_west = np.zeros(dst_v_west.shape[1])
    
    for i in range(dst_u_north.shape[1]):
        dst_ubar_north[i] = (dst_u_north[:,i] * np.diff(z_u_north[:,i])).sum() / -z_u_north[0,i]
    for i in range(dst_v_north.shape[1]):
        dst_vbar_north[i] = (dst_v_north[:,i] * np.diff(z_v_north[:,i])).sum() / -z_v_north[0,i]
    for i in range(dst_u_south.shape[1]):
        dst_ubar_south[i] = (dst_u_south[:,i] * np.diff(z_u_south[:,i])).sum() / -z_u_south[0,i]
    for i in range(dst_v_south.shape[1]):
        dst_vbar_south[i] = (dst_v_south[:,i] * np.diff(z_v_south[:,i])).sum() / -z_v_south[0,i]
    for j in range(dst_u_east.shape[1]):
        dst_ubar_east[j] = (dst_u_east[:,j] * np.diff(z_u_east[:,j])).sum() / -z_u_east[0,j]
        dst_ubar_west[j] = (dst_u_west[:,j] * np.diff(z_u_west[:,j])).sum() / -z_u_west[0,j]
    for j in range(dst_v_east.shape[1]):
        dst_vbar_east[j] = (dst_v_east[:,j] * np.diff(z_v_east[:,j])).sum() / -z_v_east[0,j]
        dst_vbar_west[j] = (dst_v_west[:,j] * np.diff(z_v_west[:,j])).sum() / -z_v_west[0,j]
        
    #mask
    dst_ubar_north = np.ma.masked_where(dst_grd.hgrid.mask_u[-1,:] == 0, dst_ubar_north)
    dst_ubar_south = np.ma.masked_where(dst_grd.hgrid.mask_u[0,:] == 0, dst_ubar_south)
    dst_ubar_east = np.ma.masked_where(dst_grd.hgrid.mask_u[:,-1] == 0, dst_ubar_east)
    dst_ubar_west = np.ma.masked_where(dst_grd.hgrid.mask_u[:,0] == 0, dst_ubar_west)
    dst_vbar_north = np.ma.masked_where(dst_grd.hgrid.mask_v[-1,:] == 0, dst_vbar_north)
    dst_vbar_south = np.ma.masked_where(dst_grd.hgrid.mask_v[0,:] == 0, dst_vbar_south)
    dst_vbar_east = np.ma.masked_where(dst_grd.hgrid.mask_v[:,-1] == 0, dst_vbar_east)
    dst_vbar_west = np.ma.masked_where(dst_grd.hgrid.mask_v[:,0] == 0, dst_vbar_west)
    
    # write data in destination file
    if verbose: print('==>Info: write data in destination file')
    ncu.variables['ocean_time'][0] = time
    ncu.variables['u_north'][0] = dst_u_north
    ncu.variables['u_south'][0] = dst_u_south
    ncu.variables['u_east'][0] = dst_u_east
    ncu.variables['u_west'][0] = dst_u_west
    ncu.variables['ubar_north'][0] = dst_ubar_north
    ncu.variables['ubar_south'][0] = dst_ubar_south
    ncu.variables['ubar_east'][0] = dst_ubar_east
    ncu.variables['ubar_west'][0] = dst_ubar_west

    ncv.variables['ocean_time'][0] = time
    ncv.variables['v_north'][0] = dst_v_north
    ncv.variables['v_south'][0] = dst_v_south
    ncv.variables['v_east'][0] = dst_v_east
    ncv.variables['v_west'][0] = dst_v_west
    ncv.variables['vbar_north'][0] = dst_vbar_north
    ncv.variables['vbar_south'][0] = dst_vbar_south
    ncv.variables['vbar_east'][0] = dst_vbar_east
    ncv.variables['vbar_west'][0] = dst_vbar_west    

    # close destination file
    src_ncuv.close(); ncu.close(); ncv.close()

