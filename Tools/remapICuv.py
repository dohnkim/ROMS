import os
import numpy as np
from netCDF4 import Dataset, date2num, num2date

import pyroms
import pyroms_toolbox

class nctime(object):
    pass

def remapICuv(src_fileuv, src_vnameuv, wts_file, src_grd, dst_grd, dst_fileuv, dxy=20, cdepth=0, kk=0, verbose=True):

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
    #Mp, Lp = dst_grd.hgrid.mask_rho.shape

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
    print('==>Info: Creating variable u')
    ncu.createVariable('u', 'f8', ('ocean_time', 's_rho', 'eta_u', 'xi_u'), fill_value=spval)
    ncu.variables['u'].long_name = '3D u-momentum component'
    ncu.variables['u'].units = 'meter second-1'
    ncu.variables['u'].field = 'u-velocity, scalar, series'
    # create variable in destination file
    print('==>Info: Creating variable ubar')
    ncu.createVariable('ubar', 'f8', ('ocean_time', 'eta_u', 'xi_u'), fill_value=spval)
    ncu.variables['ubar'].long_name = '2D u-momentum component'
    ncu.variables['ubar'].units = 'meter second-1'
    ncu.variables['ubar'].field = 'ubar-velocity,, scalar, series'

    print('==>Info: Creating variable v')
    ncv.createVariable('v', 'f8', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=spval)
    ncv.variables['v'].long_name = '3D v-momentum component'
    ncv.variables['v'].units = 'meter second-1'
    ncv.variables['v'].field = 'v-velocity, scalar, series'
    print('==>Info: Creating variable vbar')
    ncv.createVariable('vbar', 'f8', ('ocean_time', 'eta_v', 'xi_v'), fill_value=spval)
    ncv.variables['vbar'].long_name = '2D v-momentum component'
    ncv.variables['vbar'].units = 'meter second-1'
    ncv.variables['vbar'].field = 'vbar-velocity,, scalar, series'


    # remaping
    print('==>Info: remapping and rotating u and v for I.C.', ', time =', time)


    # flood the grid
    print('==>Info: flood the grid')
    src_uz = pyroms_toolbox.Grid_HYCOM.flood_fast(src_varu, src_grd, pos='t', \
                spval=spval, dxy=dxy, cdepth=cdepth, kk=kk, verbose=verbose)
    src_vz = pyroms_toolbox.Grid_HYCOM.flood_fast(src_varv, src_grd, pos='t', \
                spval=spval, dxy=dxy, cdepth=cdepth, kk=kk, verbose=verbose)

    # horizontal interpolation using scrip weights
    print('==>Info: horizontal interpolation using scrip weights')
    dst_uz = pyroms.remapping.remap(src_uz, wts_file, spval=spval)
    dst_vz = pyroms.remapping.remap(src_vz, wts_file, spval=spval)

    # vertical interpolation from standard z level to sigma
    print('==>Info: vertical interpolation from standard z level to sigma')
    dst_u = pyroms.remapping.z2roms(dst_uz[::-1,:,:], dst_grdz, \
                                    dst_grd, Cpos='rho', spval=spval, flood=False, \
                                    dmax=dxy, cdepth=cdepth, kk=kk)
    dst_v = pyroms.remapping.z2roms(dst_vz[::-1,:,:], dst_grdz,  \
                                    dst_grd, Cpos='rho', spval=spval, flood=False, \
                                    dmax=dxy, cdepth=cdepth, kk=kk)


    # rotate u,v fields
    src_angle = pyroms.remapping.remap(src_grd.angle, wts_file, spval=spval)

    dst_angle = dst_grd.hgrid.angle_rho
    angle = dst_angle - src_angle
    angle = np.tile(angle, (dst_grd.vgrid.N, 1, 1))
    U = dst_u + dst_v*1j
    eitheta = np.exp(-1j*angle[:,:,:])
    U = U * eitheta
    dst_u = np.real(U)
    dst_v = np.imag(U)


    # move back to u,v points
    dst_u = 0.5 * (dst_u[:,:,:-1] + dst_u[:,:,1:])
    dst_v = 0.5 * (dst_v[:,:-1,:] + dst_v[:,1:,:])

    # spval
    idxu = np.where(dst_grd.hgrid.mask_u == 0)
    idxv = np.where(dst_grd.hgrid.mask_v == 0)
    for n in range(dst_grd.vgrid.N):
        dst_u[n,idxu[0], idxu[1]] = spval
        dst_v[n,idxv[0], idxv[1]] = spval


    # compute depth average velocity ubar and vbar
    # get z at the right position
    z_u = 0.5 * (dst_grd.vgrid.z_w[0,:,:,:-1] + dst_grd.vgrid.z_w[0,:,:,1:])
    z_v = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,:] + dst_grd.vgrid.z_w[0,:,1:,:])

    dst_ubar = np.zeros((dst_u.shape[1], dst_u.shape[2]))
    dst_vbar = np.zeros((dst_v.shape[1], dst_v.shape[2]))

    for i in range(dst_ubar.shape[1]):
        for j in range(dst_ubar.shape[0]):
            dst_ubar[j,i] = (dst_u[:,j,i] * np.diff(z_u[:,j,i])).sum() / -z_u[0,j,i]

    for i in range(dst_vbar.shape[1]):
        for j in range(dst_vbar.shape[0]):
            dst_vbar[j,i] = (dst_v[:,j,i] * np.diff(z_v[:,j,i])).sum() / -z_v[0,j,i]

    # spval
    dst_ubar[idxu[0], idxu[1]] = spval
    dst_vbar[idxv[0], idxv[1]] = spval

    # write data in destination file
    print('==>Info: write data in destination file')
    ncu.variables['ocean_time'][0] = time
    ncu.variables['u'][0] = dst_u
    ncu.variables['ubar'][0] = dst_ubar

    ncv.variables['ocean_time'][0] = time
    ncv.variables['v'][0] = dst_v
    ncv.variables['vbar'][0] = dst_vbar

    # close destination file
    src_ncuv.close(); ncu.close(); ncv.close()

