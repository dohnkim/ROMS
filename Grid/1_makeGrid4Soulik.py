#!/usr/bin/env python
# coding: utf-8

# # Pyroms로 태풍 솔릭(Soulik) 모의를 위한 ROMS 격자 만들기
# 
# **Note**) pyroms는 육지와 바다의 masking을 마우스로 interactive하게 편집할 수 있도록 widget을 제공합니다.<br/>
# > jupyter notebook 이 아닌 jupyter lab을 사용하시는 분들은<br/>
# > "%matplotlib notebook" 이 호환이 되지 않으므로 "%matplotlib widget"를 사용하여야 합니다.<br/>
# > 이를 사용하기 위해서는 "jupyter-matplotlib"를 설치하셔야 하며 자세한 사항은  [여기](https://github.com/matplotlib/jupyter-matplotlib) 를 참조하세요.

# In[1]:


#=== for jupyter notebook ===
#%matplotlib notebook
#=== for jupyter lab ===
get_ipython().run_line_magic('matplotlib', 'widget')
#=== for general use ===
# %matplotlib inline

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata

import matplotlib.colors as colors
# from scipy.signal import medfilt2d
import netCDF4

import pyroms
import pyroms_toolbox
from bathy_smoother import *


# ## 1. 격자 영역 지정
# 
# 지도 투영 방법에 관한 개념설명은 저의 [홈페이지](http://www.dhkim.info)에 있는 "[WRF ARW에서 사용되는 지도 투영법](https://dhkim.tistory.com/273)"을 참고하세요.<br/>
# 중위도 지역의 모의에는 보통 "Lambert Conformal" 투영법을 사용하며,<br/> 
# 저위도나 적도지역은 "Mercator", 극 지역은 "Polar Stereographic",<br/>
# 전지구 모의에는 "Cylindrical Equidistant" 투영법을 사용합니다.

# In[2]:


# 격자 해상도 지정 (Lm: 경도, Mm: 위도)
Lm = 430
Mm = 320

# 위경도 영역 지정
lon0 = 117.; lat0 = 52. # 영역의 왼쪽 위
lon1 = 117.; lat1 = 20. # 영역의 왼쪽 아래
lon2 = 160.; lat2 = 20. # 영역의 오른쪽 아래
lon3 = 160.; lat3 = 52. # 영역의 오른쪽 위
centerLonLat     = [126.564, 33.457] # 영역의 중심 경위도 (투영법에서 참조하는 중심 위치임)
trueLatLowerUpper= [ 30., 40.] # 지도투영에 의해 발생하는 위도 왜곡에서 지표면과 교차되는 아랫지점과 윗지점

M = Basemap(projection='lcc', lat_0=centerLonLat[1], lon_0=centerLonLat[0],             lat_1=trueLatLowerUpper[0], lat_2=trueLatLowerUpper[1],             resolution='i',width=7000000,height=5500000) # set width & height widely just to check
            #llcrnrlon=lon1,llcrnrlat=lat1,urcrnrlon=lon3,urcrnrlat=lat3, resolution='i')
# M.latmin; M.lonmin; M.latmax; M.lonmax

# %matplotlib inline
#from matplotlib.patches import Polygon
_=M.drawcoastlines(); _=M.drawcountries(); _=M.drawmapboundary()
_=M.fillcontinents(color='coral',lake_color='aqua')
#x, y = M([lon0,lon1,lon2,lon3],[lat0,lat1,lat2,lat3])
#xy = list(zip(x,y))
#poly = Polygon(xy, facecolor='red', alpha=0.4)
#_=plt.gca().add_patch(poly)


# ## 2. pyroms로 기본 수평 격자 만들기
# 
# ### 2-1. 기본 수평 격자 만들기
# Gridgen() 함수에 인자로 들어가는 "beta" 값은 "[Pyroms - Python for ROMS](https://www.myroms.org/wiki/images/7/7f/Intro_pyroms.pdf)" 를 참고하세요.<br/>
# 영역의 초기 위치에서부터 모서리가 시계반대방향으로 꺽이면 +1 이고 시계방향으로 꺽이면 -1 입니다.

# In[3]:


lonp = np.array([lon0,lon1,lon2,lon3])
latp = np.array([lat0,lat1,lat2,lat3])
beta = np.array([1, 1, 1, 1])

## 격자의 모서리를 수정하지 않고 주어진 값으로 수평 격자를 만들 경우:
hgrd = pyroms.grid.Gridgen(lonp, latp, beta, (Mm+3, Lm+3), proj=M)

## 격자의 모서리를 직접 수정하여 수평 격자를 만들 경우:
#M.drawcoastlines()
#xp, yp = M(lonp, latp)
#bry = pyroms.hgrid.BoundaryInteractor(xp, yp, beta, shp=(Mm+3,Lm+3), proj=M)
#hgrd=bry.grd


# ### 2-2. C-grid 격자 만들기
# Gridgen()으로 만든 격자를 이용하여<br/>
# 위경도 격자체계의 Curvilinear Arakawa C-grid 를 만듭니다.

# In[4]:


lonv, latv = list(M(hgrd.x_vert, hgrd.y_vert, inverse=True))
hgrd = pyroms.grid.CGrid_geo(lonv, latv, M)

x,y = M(hgrd.lon_rho.flatten(),hgrd.lat_rho.flatten())
_=M.plot(x,y,'bo',markersize=0.1)


# ### 2-3. coastline 정보를 이용하여 육지 영역 masking 하기

# In[5]:


# generate the mask
#for verts in map.coastsegs:
#    hgrd.mask_polygon(verts)
# alternate version from johan.navarro.padron

for xx,yy in M.coastpolygons:
    xa = np.array(xx, np.float32)
    ya = np.array(yy,np.float32)
    vv = np.zeros((xa.shape[0],2))
    vv[:, 0] = xa
    vv[:, 1] = ya
    hgrd.mask_polygon(vv,mask_value=0)


# ### 2-4. masking된 육지 영역을 확인하고 수정하기
# 
# 영역을 확대 및 이동하여 masking 영역을 확인하고, 키보드의 "e"키와 마우스를 눌러서 육지와 바다 영역을 수정합니다.<br/>
# 이 기능이 정상적으로 작동하지 않는 경우는 프로그램 최상단에서 선언한<br/>
# "%matplotlib widget" 또는 "%matplotlib notebook"이 정상적으로 작동하지 않는 경우 입니다.

# In[6]:


# Edit the land mask interactively.
#pyroms.grid.edit_mask_mesh(hgrd, proj=M)
#edit_mask_mesh_ij is a faster version using imshow... but no map projection.
coast = pyroms.utility.get_coast_from_map(M)
pyroms.grid.edit_mask_mesh_ij(hgrd, coast=coast)


# ## 3. 격자에 맞는 수심자료의 내삽
# 만들어진 격자에 Etopo2 지형 자료를 사용하여 내삽을 합니다.<br/>
# Etopo2 자료는 다양한 방법으로 얻을 수있으며,<br/> 
# 여기에서는 matplotlib의 basemap에서 제공하는 자료를 사용합니다.<br/>
# 다음의 주소에서 필요한 자료를 받아서 "~/Data/Topog/ETOPO/" 에 저장하세요.<br/>
# "https://github.com/matplotlib/basemap/tree/master/examples " 

# ### 3-1. Etopo2 지형 자료 읽기
# - Etopo2 자료에서 위경도를 읽고<br/>
# - 지형자료를 읽어서 육지와 바다의 부호를 반대로 바꾸고,<br/>
# - 최소 수심을 5m 로 제한하여 육지부분과 5m 이하의 수심을 모두 5m로 제한합니다. 

# In[7]:


datadir = os.getenv('HOME')+'/Data/Topog/ETOPO/'

lats = np.loadtxt(os.path.join(datadir, 'etopo20lats.gz'))
lons = np.loadtxt(os.path.join(datadir, 'etopo20lons.gz'))

topo = np.loadtxt(os.path.join(datadir, 'etopo20data.gz'))
topo = -topo # depth positive
hmin = 5     # fix minimum depth
topo = pyroms_toolbox.change(topo, '<', hmin, hmin)


# ### 3-2. 수심 자료를 격자점에 내삽
# - 수심 자료를 격자점에 맞게 내삽하고,
# - 내삽한 수심이 최소 수심 5m를 넘지 않게 다시 한번 검사하여 수정하며,
# - 육지로 masking되어 있는 곳도 모두 5m 로 수정합니다.

# In[8]:


# interpolate new bathymetry
lon, lat = np.meshgrid(lons, lats)
#h = griddata(lon.flat,lat.flat,topo.flat,hgrd.lon_rho,hgrd.lat_rho)                     # from matplotlib
points = np.array([lon.flatten(),lat.flatten()]).T
h = griddata(np.array(points),topo.flatten(),(hgrd.lon_rho,hgrd.lat_rho),method='linear') # from scipy

# insure that depth is always deeper than hmin
h = pyroms_toolbox.change(h, '<', hmin, hmin)

# set depth to hmin where masked
idx = np.where(hgrd.mask_rho == 0)
h[idx] = hmin

# save raw bathymetry
hraw = h.copy()


# ### 3-3. 수심 평활화
# ROMS와 같이 연직 좌표계를 sigma coordinate로 사용하는 모델들은 격자 간의 수심 경사가 너무 커지지 않게 해야 합니다.<br/>
# 이를 위해서 보통 bathymetry roughness <= 0.35 정도로 제한하여 수심을 평활화 합니다.

# In[9]:


# check bathymetry roughness
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())

# smooth the raw bathy using the direct iterative method from Martinho and Batteen (2006)
rx0_max = 0.35
h = bathy_smoothing.smoothing_Positive_rx0(hgrd.mask_rho, h, rx0_max)

# check bathymetry roughness again
RoughMat = bathy_tools.RoughnessMatrix(h, hgrd.mask_rho)
print('Max Roughness value is: ', RoughMat.max())


# ## 4. 연직좌표계 결정 및 Grid file 저장

# In[10]:


# vertical coordinate
theta_b = 0.1
theta_s = 7.0
Tcline = 50
N = 30
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)


# In[11]:


get_ipython().run_line_magic('env', 'PYROMS_GRIDID_FILE=./gridid.txt')

grd_name = 'Soulik'; version = 'v2'
grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

# write grid to netcdf file
pyroms.grid.write_ROMS_grid(grd, filename='Soulik_grd_'+version+'.nc')

