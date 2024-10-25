
import netCDF4
import numpy as np
from operator import add
from datetime import datetime
from arpifs4py import wlfa as LFA4py
import os
from functools import reduce

os.putenv('DR_HOOK_NOT_MPI','1')  # <== Necessaire pour ne pas partir dans des routines MPI

#******************************************
#  Definition de la classe NetCdf         #
#******************************************
class NetCdf:


#-----------------------------------------------------------------------------------#        
   def __init__(self,full=False):

     
# some usefull constants ...   
     self.Rg      = 9.80665
     self.Rtt     = 273.16
     self.Cpd     = 1004.7
     self.Cpv     = 1846.1
     self.Rd      = 287.06
     self.Rv      = 461.52
     self.Rcw     = 4218.
     self.Rcs     = 2106.
     self.Rlvtt   = 2.5008E+6
     self.Rlstt   = 2.8345E+6
     self.Rlvzer  = self.Rlvtt + self.Rtt*(self.Rcw-self.Cpv)
     self.Rlszer  = self.Rlstt + self.Rtt*(self.Rcs-self.Cpv)
     self.P0      = 101325.
     self.sigma   = 5.670509E-8
     self.epsilon = 0.965
     self.RA      = 6378.137   
     
     self.full    = full

#    Variables intensives 
     self.vuui_list = []
     self.vvvi_list = []
     self.vqvi_list = []
     self.vqli_list = []
     self.vhri_list = []
     self.vqni_list = []
     self.vnti_list = []
     self.vcti_list = []
     self.vomi_list = []
     self.vep_list  = []
     self.DeltaP_sur_g_list = []
     self.DeltaP_list = []
        
#    Initialisation des champs de sortie
     self.nlev_list   = []
     self.nlevp1_list = []
     self.P_half_list = []
     self.P_full_list = []
     self.vuu_list = []
     self.vvv_list = []
     self.dd_list  = []
     self.ff_list  = []
     self.vqv_list = []
     self.vql_list = []
     self.vhr_list = []
     
     self.vqn_list = []
     self.vnt_list = []
     self.vzf_list = []
     self.vzh_list = []
     self.vct_list = []
     self.theta_list = []
     self.thetav_list = []
     self.vom_list = []
     self.Phi_surf_list = []
     self.time_list = [0]
     self.timef_list = []
     self.Ps_list   = []
     self.albedo_list = []
     self.gz0d_list = []
     self.gz0t_list = []
     self.pblh_list = []
     self.tsurf_list = [0]
     self.t2m_list = [0]
     self.q2m_list = [0]
     self.hu2m_list = [0]
     self.u10m_list = [0]
     self.v10m_list = [0]
     self.SWdC_list = []
     self.LWdC_list = []
     self.SWnC_list = []
     self.LWnC_list = []
     self.LheatC_list = []
     self.SheatC_list = []
     self.FpreciplC_list = []
     self.FprecipsC_list = []
     self.FpreciplconC_list = []
     self.FprecipsconC_list = []
     self.FpreciplstrC_list = []
     self.FprecipsstrC_list = []
     self.FraySolC_list = []
     self.FrayTerC_list = []
     self.Pacc_list  = []
     self.Psnowacc_list = []
     self.Snow_list = []
     self.Pl_list = []
     self.Psnow_list = []
     self.Ustressacc_list = []
     self.Vstressacc_list = []
     self.lsm_list = []
     self.t_soil_h = []
     self.tsrad_h = []
# Advection
     self.tctresiC_list = []    
     self.tqvresiC_list = []    
     self.tuuresiC_list = []    
     self.tvvresiC_list = []     
                
     self.l_tctresi = 0
     self.l_tqvresi = 0
     self.l_tuuresi = 0
     self.l_tvvresi = 0     

# Qv tendencies => cumulated flux, as red in ddh files

     self.qfluxC_trans_turb_list  = []   # Turbulence transport tendency <= Ce sont des flux en fait, transformé en tendance ailleurs
     self.qfluxC_trans_conv_list  = []   # Convective transport tendency
     self.qfluxC_cond_conv_l_list = []   # Convective condensation liquid
     self.qfluxC_cond_conv_i_list = []   # Convective condensation ice
     self.qfluxC_cond_stra_l_list = []   # Stratiforme condensation liquid
     self.qfluxC_cond_stra_i_list = []   # Stratiforme condensation ice
     self.qfluxC_evap_stra_l_list = []   # Stratiforme evaporation liquid
     self.qfluxC_evap_stra_i_list = []   # Stratiforme evaporation ice

#  H tendencies 

     self.hfluxC_solar_list         = []   # Solar H flux
     self.hfluxC_therm_list         = []   # Thermal H flux
     self.hfluxC_turb_list          = []   # Turbulent H flux
     self.hfluxC_conv_list          = []   # Convective H flux

#  U and V tendencies

     self.UfluxC_conv_list          = []     
     self.VfluxC_conv_list          = []     
     self.UfluxC_turb_list          = []     
     self.VfluxC_turb_list          = []     
     self.UfluxC_gwd_list           = []     
     self.VfluxC_gwd_list           = []     
     
#-----------------------------------------------------------------------------------#        
   def __add__(self,add_object):

     return self.oper(add_object)

#-----------------------------------------------------------------------------------#        
   def __truediv__(self,scalar):

     return self.oper(scalar)

#-----------------------------------------------------------------------------------#        
   def oper(self,oper):

# Aditionne 2 objets ou divise un objet par un float
     out_fields = NetCdf()

     out_fields.nlev       = self.nlev
     out_fields.nlev_list  = self.nlev_list
     out_fields.nlevp1     = self.nlevp1
     out_fields.nlevp1_list= self.nlevp1_list
     out_fields.time_slot  = self.time_slot
     out_fields.time_list  = self.time_list
     out_fields.timef_list = self.timef_list
     out_fields.time_unit  = self.time_unit
     out_fields.nsurf      = self.nsurf
     out_fields.model      = self.model
     out_fields.point_name = 'point_name'
     out_fields.file_date  = self.file_date
     out_fields.version    = self.version
     out_fields.A_list     = self.A_list
     out_fields.B_list     = self.B_list
     out_fields.ddt        = self.ddt
     out_fields.hour       = self.hour

#  On moyenne les positions des points

     if type(oper) == float :
        scalar = oper
        out_fields.lat        = self.lat  / scalar
        out_fields.lon        = self.lon  / scalar
        out_fields.latr       = self.latr / scalar
        out_fields.lonr       = self.lonr / scalar
        out_fields.dist       = self.dist / scalar
     else:   
        add_object = oper
        out_fields.lat        = self.lat  + add_object.lat
        out_fields.lon        = self.lon  + add_object.lon
        out_fields.latr       = self.latr + add_object.latr
        out_fields.lonr       = self.lonr + add_object.lonr
        out_fields.dist       = self.dist + add_object.dist

     if self.nsurf != 0:
        out_fields.SWd_list = []  
        out_fields.LWd_list = [] 
        out_fields.SWn_list = []  
        out_fields.LWn_list = []  
        out_fields.Lheat_list = []  
        out_fields.Sheat_list = [] 
        out_fields.t_soil_h = [] 
        out_fields.tsrad_h = []
   #    if self.nsurf == 1:
        out_fields.Fprecipl_list = [] 
        out_fields.Fprecips_list = [] 
        out_fields.Ustress_list = [] 
        out_fields.Vstress_list = [] 

# Division par un scalaire ou addition d'un autre objet
     time = 0
     while time <= self.time_slot:
        vctloc_list   = []
        vqvloc_list   = []
        vqlloc_list   = []
        vhrloc_list   = []
        vqnloc_list   = []
        vntloc_list   = []
        veploc_list   = []
        vomloc_list   = []
        vdpsgloc_list = []
        vuuloc_list   = []
        vvvloc_list   = []

# Surface fields
        if self.nsurf != 0:
          if type(oper) == float :
             scalar = oper
             Phi_surfloc_list = self.Phi_surf_list[time] / scalar
             lsmloc_list      = self.lsm_list[time]      / scalar
             Snowloc_list     = self.Snow_list[time]     / scalar
             t_soil_hloc      = self.t_soil_h[time]      / scalar
             tsrad_hloc      = self.tsrad_h[time]      / scalar
          else:   
             add_object = oper
             Phi_surfloc_list = self.Phi_surf_list[time] + add_object.Phi_surf_list[time]
             lsmloc_list      = self.lsm_list[time]      + add_object.lsm_list[time]
             Snowloc_list     = self.Snow_list[time]     + add_object.Snow_list[time]
             t_soil_hloc      = self.t_soil_h[time]      + add_object.t_soil_h[time]
             tsrad_hloc      = self.tsrad_h[time]      + add_object.tsrad_h[time]
             
        lev           = 0
        while lev < self.nlev:
           if type(oper) == float :
             scalar = oper
             vctloc_list.append(self.vcti_list[time][lev] / scalar)
             vqvloc_list.append(self.vqvi_list[time][lev] / scalar)
             vqlloc_list.append(self.vqli_list[time][lev] / scalar)
             vhrloc_list.append(self.vhri_list[time][lev] / scalar)

             vqnloc_list.append(self.vqni_list[time][lev] / scalar)
             vntloc_list.append(self.vnti_list[time][lev] / scalar)
             veploc_list.append(self.vep_list[time][lev]  / scalar)
             vomloc_list.append(self.vomi_list[time][lev] / scalar)
             vuuloc_list.append(self.vuui_list[time][lev] / scalar)
             vvvloc_list.append(self.vvvi_list[time][lev] / scalar)
             vdpsgloc_list.append(self.DeltaP_sur_g_list[time][lev] / scalar)
             

           else :  
             add_object = oper
             vctloc_list.append(self.vcti_list[time][lev] + add_object.vcti_list[time][lev])
             vqvloc_list.append(self.vqvi_list[time][lev] + add_object.vqvi_list[time][lev])
             vqlloc_list.append(self.vqli_list[time][lev] + add_object.vqli_list[time][lev])
             vhrloc_list.append(self.vhri_list[time][lev] + add_object.vhri_list[time][lev])

             vqnloc_list.append(self.vqni_list[time][lev] + add_object.vqni_list[time][lev])
             vntloc_list.append(self.vnti_list[time][lev] + add_object.vnti_list[time][lev])
             veploc_list.append(self.vep_list[time][lev]  + add_object.vep_list[time][lev])
             vomloc_list.append(self.vomi_list[time][lev] + add_object.vomi_list[time][lev])
             vuuloc_list.append(self.vuui_list[time][lev] + add_object.vuui_list[time][lev])
             vvvloc_list.append(self.vvvi_list[time][lev] + add_object.vvvi_list[time][lev])
             vdpsgloc_list.append(self.DeltaP_sur_g_list[time][lev] + add_object.DeltaP_sur_g_list[time][lev])
             
           lev += 1           
           
        out_fields.vcti_list.append(vctloc_list)
        out_fields.vqvi_list.append(vqvloc_list)
        out_fields.vqli_list.append(vqlloc_list)
        out_fields.vhri_list.append(vhrloc_list)

        out_fields.vqni_list.append(vqnloc_list)
        out_fields.vnti_list.append(vntloc_list)
        out_fields.vep_list.append(veploc_list)
        out_fields.vomi_list.append(vomloc_list)
        out_fields.vuui_list.append(vuuloc_list)
        out_fields.vvvi_list.append(vvvloc_list)
        out_fields.DeltaP_sur_g_list.append(vdpsgloc_list)
        
        if self.nsurf != 0:
           out_fields.Phi_surf_list.append(Phi_surfloc_list)
           out_fields.lsm_list.append(lsmloc_list)
           out_fields.Snow_list.append(Snowloc_list)
           out_fields.t_soil_h.append(t_soil_hloc)
           out_fields.tsrad_h.append(tsrad_hloc)

        time += 1    


# Champs de surface aux temps intermediaires ...
     if self.nsurf != 0:
        time = 0
        while time < self.time_slot:

           if type(oper) == float :
              scalar = oper
              gz0dloc_list     = self.gz0d_list[time]     / scalar
              gz0tloc_list     = self.gz0t_list[time]     / scalar
              pblhloc_list     = self.pblh_list[time]     / scalar
              SWdloc_list      = self.SWd_list[time]      / scalar
              LWdloc_list      = self.LWd_list[time]      / scalar
              albedoloc_list   = self.albedo_list[time]   / scalar
              SWnloc_list      = self.SWn_list[time]      / scalar
              LWnloc_list      = self.LWn_list[time]      / scalar
              Lheatloc_list    = self.Lheat_list[time]    / scalar
              Sheatloc_list    = self.Sheat_list[time]    / scalar
              Paccloc_list     = self.Pacc_list[time]     / scalar
              Psnowaccloc_list = self.Psnowacc_list[time] / scalar
              tsurfloc_list      = self.tsurf_list[time]      / scalar
              t2mloc_list      = self.t2m_list[time]      / scalar
              q2mloc_list      = self.q2m_list[time]      / scalar
              if self.nsurf == 2: hu2mloc_list     = self.hu2m_list[time]      / scalar
              u10mloc_list     = self.u10m_list[time]     / scalar
              v10mloc_list     = self.v10m_list[time]     / scalar
              Plloc_list       = self.Pl_list[time]       / scalar
              Psnowloc_list    = self.Psnow_list[time]    / scalar
    #         if self.nsurf == 1:
              Fpreciplloc_list  = self.Fprecipl_list[time]  / scalar
              Fprecipsloc_list  = self.Fprecips_list[time]  / scalar
              Ustressloc_list   = self.Ustress_list[time]   / scalar
              Vstressloc_list   = self.Vstress_list[time]   / scalar

           else :  
              add_object = oper
              gz0dloc_list     = self.gz0d_list[time]     + add_object.gz0d_list[time]
              gz0tloc_list     = self.gz0t_list[time]     + add_object.gz0t_list[time]
              pblhloc_list     = self.pblh_list[time]     + add_object.pblh_list[time]
              tsradloc_list     = self.tsrad_list[time]     + add_object.tsrad_list[time]
              SWdloc_list      = self.SWd_list[time]      + add_object.SWd_list[time]
              LWdloc_list      = self.LWd_list[time]      + add_object.LWd_list[time]
              albedoloc_list   = self.albedo_list[time]   + add_object.albedo_list[time]
              SWnloc_list      = self.SWn_list[time]      + add_object.SWn_list[time]
              LWnloc_list      = self.LWn_list[time]      + add_object.LWn_list[time]
              Lheatloc_list    = self.Lheat_list[time]    + add_object.Lheat_list[time]
              Sheatloc_list    = self.Sheat_list[time]    + add_object.Sheat_list[time]
              Paccloc_list     = self.Pacc_list[time]     + add_object.Pacc_list[time]
              Psnowaccloc_list = self.Psnowacc_list[time] + add_object.Psnowacc_list[time]
              tsurfloc_list      = self.tsurf_list[time]      + add_object.tsurf_list[time]
              t2mloc_list      = self.t2m_list[time]      + add_object.t2m_list[time]
              q2mloc_list      = self.q2m_list[time]      + add_object.q2m_list[time]
              if self.nsurf == 2: hu2mloc_list     = self.hu2m_list[time]     + add_object.hu2m_list[time]
              u10mloc_list     = self.u10m_list[time]     + add_object.u10m_list[time]
              v10mloc_list     = self.v10m_list[time]     + add_object.v10m_list[time]
              Plloc_list       = self.Pl_list[time]       + add_object.Pl_list[time]
              Psnowloc_list    = self.Psnow_list[time]    + add_object.Psnow_list[time]
             
    #         if self.nsurf == 1:
              Fpreciplloc_list  = self.Fprecipl_list[time]  + add_object.Fprecipl_list[time]
              Fprecipsloc_list  = self.Fprecips_list[time]  + add_object.Fprecips_list[time]
              Ustressloc_list   = self.Ustress_list[time]   + add_object.Ustress_list[time]
              Vstressloc_list   = self.Vstress_list[time]   + add_object.Vstress_list[time]
             

           out_fields.gz0d_list.append(gz0dloc_list)
           out_fields.gz0t_list.append(gz0tloc_list)
           out_fields.pblh_list.append(pblhloc_list)
           out_fields.tsrad_list.append(tsradloc_list)
           out_fields.SWd_list.append(SWdloc_list)
           out_fields.LWd_list.append(LWdloc_list)
           out_fields.albedo_list.append(albedoloc_list)
           out_fields.SWn_list.append(SWnloc_list)
           out_fields.LWn_list.append(LWnloc_list)
           out_fields.Lheat_list.append(Lheatloc_list)
           out_fields.Sheat_list.append(Sheatloc_list)
           out_fields.Pacc_list.append(Paccloc_list)
           out_fields.Psnowacc_list.append(Psnowaccloc_list)
           out_fields.tsurf_list.append(tsurfloc_list)
           out_fields.t2m_list.append(t2mloc_list)
           out_fields.q2m_list.append(q2mloc_list)
           if self.nsurf == 2: out_fields.hu2m_list.append(hu2mloc_list)
           out_fields.u10m_list.append(u10mloc_list)
           out_fields.v10m_list.append(v10mloc_list)
           out_fields.Pl_list.append(Plloc_list)
           out_fields.Psnow_list.append(Psnowloc_list)
           
   #       if self.nsurf == 1:
           out_fields.Fprecipl_list.append(Fpreciplloc_list)
           out_fields.Fprecips_list.append(Fprecipsloc_list)
           out_fields.Ustress_list.append(Ustressloc_list)
           out_fields.Vstress_list.append(Vstressloc_list)

           time += 1    

     if self.nsurf != 0:
       out_fields.tsrad_list[0] = out_fields.tsrad_list[1]
       out_fields.tsurf_list[0] = out_fields.tsurf_list[1]
       out_fields.t2m_list[0] = out_fields.t2m_list[1]
       out_fields.q2m_list[0] = out_fields.q2m_list[1]
       if self.nsurf == 2: out_fields.hu2m_list[0] = out_fields.hu2m_list[1]
       out_fields.u10m_list[0] = out_fields.u10m_list[1]
       out_fields.v10m_list[0] = out_fields.v10m_list[1] 
     return out_fields

#-----------------------------------------------------------------------------------#        
   def convert(self):
     time = 0
     
     # Reinitialisation des champs de sortie

     self.vqv_list    = []  
     self.vql_list    = []  
     self.vhr_list    = []  

     self.vqn_list    = []   
     self.vct_list    = []        
     self.vnt_list    = []
     self.vzf_list    = []
     self.vzh_list    = []
     self.vom_list    = []        
     self.vuu_list    = []
     self.vvv_list    = []
     self.dd_list     = []
     self.ff_list     = []
     self.Ps_list     = []
     self.theta_list  = []
     self.thetav_list = []
     self.P_half_list = []
     self.P_full_list = []

     while time <= self.time_slot:
        vqvloc_list   = []
        vqlloc_list   = []
        vhrloc_list   = []

        vqnloc_list   = []
        vctloc_list   = []
        vuuloc_list   = []
        vvvloc_list   = []
        vntloc_list   = []
        vzfloc_list   = []
        vzhloc_list   = []
        vomloc_list   = []
        pshloc        = 0
        psfloc        = 0
        phalfloc_list = [0]
        pfullloc_list = []
        thetaloc_list = []
        thetavloc_list= []
        gDeltaZloc_list= []
        lev  = 0
        while lev < self.nlev:
           DeltaP_sur_g = self.DeltaP_sur_g_list[time][lev]
           DeltaP       = DeltaP_sur_g*self.Rg
           pshloc       = pshloc + DeltaP
           psfloc       = pshloc - DeltaP/2.
           DeltaP_sur_P = DeltaP / psfloc 
           phalfloc_list.append(pshloc)
           pfullloc_list.append(psfloc)
           vqv = self.vqvi_list[time][lev] / DeltaP_sur_g
           vqvloc_list.append(vqv)
           vqlloc_list.append(self.vqli_list[time][lev] / DeltaP_sur_g)
           vhrloc_list.append(self.vhri_list[time][lev] / DeltaP_sur_g)

           vqnloc_list.append(self.vqni_list[time][lev] / DeltaP_sur_g)
           vuuloc_list.append(self.vuui_list[time][lev] / DeltaP_sur_g)
           vvvloc_list.append(self.vvvi_list[time][lev] / DeltaP_sur_g)           
           vntloc_list.append(self.vnti_list[time][lev] / DeltaP_sur_g)
           vzfloc_list.append(self.vep_list[time][lev]  / DeltaP)
           vomloc_list.append(self.vomi_list[time][lev] / DeltaP_sur_g)
           Cp       = self.Cpd + (self.Cpv-self.Cpd)*vqv
           R        = self.Rd  + (self.Rv - self.Rd)*vqv
           R_sur_cp = R/Cp
           
           vct = self.vcti_list[time][lev] / DeltaP_sur_g / Cp
           RT  = R * vct
           vctloc_list.append(vct)
           gDeltaZloc_list.append(DeltaP_sur_P * RT)
           theta = vct*(self.P0/psfloc)**R_sur_cp
           thetaloc_list.append(theta)
           thetavloc_list.append(theta*R/self.Rd)
           lev += 1

        
        
        Z_top    = sum(gDeltaZloc_list)  / self.Rg
        Zh       = [Z_top]
        ilev=1
        while ilev < self.nlevp1:
           Zh.append(Zh[ilev-1]-gDeltaZloc_list[ilev-1]/self.Rg)
           ilev = ilev+1
        
        ddff = self.ddff_from_u_and_v(np.array(vuuloc_list),np.array(vvvloc_list))   
        self.vqv_list.append(vqvloc_list)  
        self.vql_list.append(vqlloc_list)  
        self.vhr_list.append(vhrloc_list)  

        self.vqn_list.append(vqnloc_list)   
        self.vct_list.append(vctloc_list)        
        self.vnt_list.append(vntloc_list)
        self.vzf_list.append(vzfloc_list)
        self.vzh_list.append(Zh)
        self.vom_list.append(vomloc_list)        
        self.vuu_list.append(vuuloc_list)
        self.vvv_list.append(vvvloc_list)
        self.dd_list.append(ddff[1])
        self.ff_list.append(ddff[0])
        self.Ps_list.append(pshloc)
        self.theta_list.append(thetaloc_list)
        self.thetav_list.append(thetavloc_list)
        self.P_half_list.append(phalfloc_list)
        self.P_full_list.append(pfullloc_list)
        time += 1

#-----------------------------------------------------------------------------------#        
   def filename(self):   
     self.file_name = self.model + '-L' + str(self.nlev) + '_' + self.point_name + '_' + self.file_date + '.nc'

#-----------------------------------------------------------------------------------#        

   def init(self,ddh_list,number,name,model,pl,debug,version):

# Initialize software version
     self.version = version
     self.model   = model
     self.debug   = debug
# Verification of the validity of the point name

     file_test = ddh_list[0]
     q_unit = LFA4py.wlfaouv(file_test,'R')
     docfichier = self.read_lfa_field(q_unit,'DOCFICHIER')
     self.nbpoint = docfichier[14]
     if pl: print('Nb points : ',self.nbpoint)
     self.nlev = docfichier[5]
     self.nlevp1 = self.nlev+1
     ddt = self.read_lfa_field(q_unit,'DATE')
     year        = ddt[0]
     month       = ddt[1]
     day         = ddt[2]
     self.hour   = ddt[3]
     minute      = ddt[4]
     processtype = ddt[8]
     self.ddt = year*10000 + month*100 + day
     self.file_date = str(self.ddt) + '{:02d}'.format(self.hour)     
     self.time_unit = 'seconds since ' + '{:04d}'.format(year) + '-' + '{:02d}'.format(month) + '-' + '{:02d}'.format(day) \
     + ' ' + '{:02d}'.format(self.hour) + ':' + '{:02d}'.format(minute) + ':00'     
     ipoint = 1
     self.extract = 0   
     if number == 0 or pl :      # <== Extraction par nom
       while ipoint <= self.nbpoint:
          try:
            (self.lon,self.lat,self.lonr,self.latr,self.dist) = self.position(q_unit,ipoint)
            self.point_name = self.point_name_maker(self.latr,self.lonr)
            if pl: print('point number : {:03d}'.format(ipoint),' latr = {:.1f}'.format(self.latr),' lat = {:.1f}'.format(self.lat),\
            ' lonr = {:.1f}'.format(self.lonr),' lon = {:.1f}'.format(self.lon),' dist = {:.1f}'.format(self.dist),' name :',self.point_name)
            if self.point_name == name:
               self.extract = 1
               self.filename()
               self.point = ipoint
               break
          except TypeError:
            print('Bad data type return from routine position. Documentation field seems to be missing')    
               
          ipoint = ipoint+1
     else :               # <== Extraction par numero
       try:
         (self.lon,self.lat,self.lonr,self.latr,self.dist) = self.position(q_unit,number)
       except TypeError:
         print('Bad data type return from routine position. Documentation field seems to be missing')     
         raise RuntimeError
       finally :
         self.point = number     
         self.point_name = name
         self.extract = 1
         self.filename()
             


#  test nomenclature ARPEGE ou AROME pour les champs de surface
     self.nsurf = 0   # Pas de champs de surface
     try:
       (fieldtype, fieldlength)  = LFA4py.wlfacas(q_unit, 'SVGFS01') 
       if fieldtype == 'R4': self.nsurf = 1  # Arpege case
     except RuntimeError:
       if self.debug : print('SVGFS01 not found in file')
     try:
       (fieldtype, fieldlength)  = LFA4py.wlfacas(q_unit, 'SVTCLS')
       if fieldtype == 'R4': self.nsurf = 2  # Arome case
     except RuntimeError:
       if self.debug : print('SVTCLS not found in file')

#  Test de la presence des champs d'advection     
     try:
       (fieldtype, fieldlength)  = LFA4py.wlfacas(q_unit, 'TCTRESI')
       if fieldtype == 'R4': self.l_tctresi = 1
     except RuntimeError:
       if self.debug : print('TCTRESI not found in file')
     try:
       (fieldtype, fieldlength)  = LFA4py.wlfacas(q_unit, 'TQVRESI')
       if fieldtype == 'R4': self.l_tqvresi = 1
     except RuntimeError:
       if self.debug : print('TQVRESI not found in file')
     try:  
       (fieldtype, fieldlength)  = LFA4py.wlfacas(q_unit, 'TUURESI')
       if fieldtype == 'R4': self.l_tuuresi = 1
     except RuntimeError:
       if self.debug : print('TUURESI not found in file')
     try:  
       (fieldtype, fieldlength)  = LFA4py.wlfacas(q_unit, 'TVVRESI')
       if fieldtype == 'R4': self.l_tvvresi = 1
     except RuntimeError:
       if self.debug : print('TVVRESI not found in file')

        
     LFA4py.wlfafer(q_unit)   

#-----------------------------------------------------------------------------------#        
   def flux_2_tend(self,flux):
   
 #     tend_list = []
      tend = np.zeros((self.time_slot,self.nlev))
      npflux = np.array(flux)
      npdeltap = np.array(self.DeltaP_sur_g_list)
   
      time = 0
      while time<self.time_slot :
         ilev = 0
         while ilev<self.nlev :
            dpog = 0.5*(npdeltap[time][ilev] + npdeltap[time+1][ilev])
            tend[time][ilev] = (npflux[time][ilev] - npflux[time][ilev+1])/dpog
            ilev+=1
         time+=1   
      return tend  
#-----------------------------------------------------------------------------------#        
   def var_2_tend(self,var):
   
      tend  = np.zeros((self.time_slot,self.nlev))
      npvar = np.array(var)
   
      time = 0
      while time<self.time_slot :
         ilev = 0
         while ilev<self.nlev :
            tend[time][ilev] = (npvar[time+1][ilev] - npvar[time][ilev])*24.
            ilev+=1
         time+=1   
      return tend  
#-----------------------------------------------------------------------------------#        
   def var_2_varf(self,var):
   
      varf  = np.zeros((self.time_slot,self.nlev+1))
      npvar = np.array(var)
   
      time = 0
      while time<self.time_slot :
         varf[time][0] = 0.5*(var[time][0] + var[time+1][0])
         varf[time][self.nlev] = 0.5*(var[time][self.nlev-1] + var[time+1][self.nlev-1])
         ilev = 1
         while ilev<self.nlev :
            varf[time][ilev] = 0.25*(var[time][ilev] + var[time+1][ilev] + var[time][ilev-1] + var[time+1][ilev-1])
            ilev+=1
         time+=1   
      return varf
#-----------------------------------------------------------------------------------#        
   def read_lfa_file(self,ddh_list):

        self.maxdps = 0.
        self.tmax   = 0
        self.time_slot = len(ddh_list)
        time = 0
        self.read_lfa_file_main(ddh_list,time,'0')
        while time < self.time_slot:
           self.read_lfa_file_main(ddh_list,time,'1')
           time = time+1
                      
        self.Phi_surf_list[0] = self.Phi_surf_list[1]   
        
        if self.nsurf != 0:
# On reporte la temperature a 1h a l'instant initial
          self.tsurf_list[0] = self.tsurf_list[1]
          self.t2m_list[0] = self.t2m_list[1]
          self.q2m_list[0] = self.q2m_list[1]
          if self.nsurf == 2: 
            self.hu2m_list[0] = self.hu2m_list[1]
          self.u10m_list[0] = self.u10m_list[1]
          self.v10m_list[0] = self.v10m_list[1]
        

# Decumul des flux et initialisation des niveaux 
          self.SWd_list = []  ; self.SWd_list.append(self.SWdC_list[0]/3600.)
          self.LWd_list = []  ; self.LWd_list.append(self.LWdC_list[0]/3600.)
          self.SWn_list = []  ; self.SWn_list.append(self.SWnC_list[0]/3600.)
          self.LWn_list = []  ; self.LWn_list.append(self.LWnC_list[0]/3600.)
          self.LWu_list = []  ; self.LWu_list.append(self.LWd_list[0] - self.LWn_list[0])
          self.SWu_list = []  ; self.SWu_list.append(self.SWd_list[0] - self.SWn_list[0])
          self.Lheat_list = []  ; self.Lheat_list.append(self.LheatC_list[0]/3600.)
          self.Sheat_list = []  ; self.Sheat_list.append(self.SheatC_list[0]/3600.)
    #     if self.nsurf == 1:
          self.Fprecipl_list = []  ; self.Fprecipl_list.append(self.FpreciplC_list[0]*24.)
          self.Fprecips_list = []  ; self.Fprecips_list.append(self.FprecipsC_list[0]*24.)
          if self.full :
            self.Fpreciplcon_list = []  ; self.Fpreciplcon_list.append(self.FpreciplconC_list[0]*24.)
            self.Fprecipscon_list = []  ; self.Fprecipscon_list.append(self.FprecipsconC_list[0]*24.)
            self.Fpreciplstr_list = []  ; self.Fpreciplstr_list.append(self.FpreciplstrC_list[0]*24.)
            self.Fprecipsstr_list = []  ; self.Fprecipsstr_list.append(self.FprecipsstrC_list[0]*24.)
          self.Ustress_list = [] ; self.Ustress_list.append(self.Ustressacc_list[0]/3600.)
          self.Vstress_list = [] ; self.Vstress_list.append(self.Vstressacc_list[0]/3600.)

          self.FraySol_list = [] ; self.FraySol_list.append(self.FraySolC_list[0]/3600.)
          self.FrayTer_list = [] ; self.FrayTer_list.append(self.FrayTerC_list[0]/3600.)
    
          self.Pl_list = []  ; self.Pl_list.append(self.Pacc_list[0]*24.)
          self.Psnow_list = []  ; self.Psnow_list.append(self.Psnowacc_list[0]*24.)

          if self.full :
# Qv tendencies => cumulated flux, as red in ddh files
            self.qflux_trans_turb_list = [] ; self.qflux_trans_turb_list.append(self.qfluxC_trans_turb_list[0]*24.)
            self.qflux_trans_conv_list = [] ; self.qflux_trans_conv_list.append(self.qfluxC_trans_conv_list[0]*24.)
            self.qflux_cond_conv_l_list = [] ; self.qflux_cond_conv_l_list.append(self.qfluxC_cond_conv_l_list[0]*24.)
            self.qflux_cond_conv_i_list = [] ; self.qflux_cond_conv_i_list.append(self.qfluxC_cond_conv_i_list[0]*24.)
            self.qflux_cond_stra_l_list = [] ; self.qflux_cond_stra_l_list.append(self.qfluxC_cond_stra_l_list[0]*24.)
            self.qflux_cond_stra_i_list = [] ; self.qflux_cond_stra_i_list.append(self.qfluxC_cond_stra_i_list[0]*24.)
            self.qflux_evap_stra_l_list = [] ; self.qflux_evap_stra_l_list.append(self.qfluxC_evap_stra_l_list[0]*24.)
            self.qflux_evap_stra_i_list = [] ; self.qflux_evap_stra_i_list.append(self.qfluxC_evap_stra_i_list[0]*24.)
          
# H tendencies     

            self.hflux_solar_list = [] ; self.hflux_solar_list.append(self.hfluxC_solar_list[0]*24.)
            self.hflux_therm_list = [] ; self.hflux_therm_list.append(self.hfluxC_therm_list[0]*24.)
            self.hflux_turb_list = [] ; self.hflux_turb_list.append(self.hfluxC_turb_list[0]*24.)
            self.hflux_conv_list = [] ; self.hflux_conv_list.append(self.hfluxC_conv_list[0]*24.)
     
#  U and V tendencies

            self.Uflux_conv_list = [] ; self.Uflux_conv_list.append(self.UfluxC_conv_list[0]*24.)
            self.Vflux_conv_list = [] ; self.Vflux_conv_list.append(self.VfluxC_conv_list[0]*24.)
            self.Uflux_turb_list = [] ; self.Uflux_turb_list.append(self.UfluxC_turb_list[0]*24.)
            self.Vflux_turb_list = [] ; self.Vflux_turb_list.append(self.VfluxC_turb_list[0]*24.)
            self.Uflux_gwd_list = []  ; self.Uflux_gwd_list.append(self.UfluxC_gwd_list[0]*24.)
            self.Vflux_gwd_list = []  ; self.Vflux_gwd_list.append(self.VfluxC_gwd_list[0]*24.)
          
          if self.l_tctresi: 
             self.tctresi_list = []
             self.tctresi_list.append(self.tctresiC_list[0]/3600.)
          if self.l_tqvresi: 
             self.tqvresi_list = []
             self.tqvresi_list.append(self.tqvresiC_list[0]/3600.)
          if self.l_tuuresi: 
             self.tuuresi_list = []
             self.tuuresi_list.append(self.tuuresiC_list[0]/3600.)
          if self.l_tvvresi: 
             self.tvvresi_list = []
             self.tvvresi_list.append(self.tvvresiC_list[0]/3600.)
          time = 1   
          while time < self.time_slot:  
             self.SWd_list.append((self.SWdC_list[time]-self.SWdC_list[time-1])/3600.) 
             self.LWd_list.append((self.LWdC_list[time]-self.LWdC_list[time-1])/3600.) 
             self.SWn_list.append((self.SWnC_list[time]-self.SWnC_list[time-1])/3600.) 
             self.LWn_list.append((self.LWnC_list[time]-self.LWnC_list[time-1])/3600.) 
             self.LWu_list.append(self.LWd_list[time] - self.LWn_list[time])
             self.SWu_list.append(self.SWd_list[time] - self.SWn_list[time])
             self.Lheat_list.append((self.LheatC_list[time]-self.LheatC_list[time-1])/3600.) 
             self.Sheat_list.append((self.SheatC_list[time]-self.SheatC_list[time-1])/3600.) 
      #      if self.nsurf == 1:
             self.Fprecipl_list.append((self.FpreciplC_list[time]-self.FpreciplC_list[time-1])*24.) 
             self.Fprecips_list.append((self.FprecipsC_list[time]-self.FprecipsC_list[time-1])*24.) 
             if self.full :
               self.Fpreciplcon_list.append((self.FpreciplconC_list[time]-self.FpreciplconC_list[time-1])*24.) 
               self.Fprecipscon_list.append((self.FprecipsconC_list[time]-self.FprecipsconC_list[time-1])*24.) 
               self.Fpreciplstr_list.append((self.FpreciplstrC_list[time]-self.FpreciplstrC_list[time-1])*24.) 
               self.Fprecipsstr_list.append((self.FprecipsstrC_list[time]-self.FprecipsstrC_list[time-1])*24.) 
             self.Ustress_list.append((self.Ustressacc_list[time]-self.Ustressacc_list[time-1])/3600.)
             self.Vstress_list.append((self.Vstressacc_list[time]-self.Vstressacc_list[time-1])/3600.)

             self.FraySol_list.append((self.FraySolC_list[time]-self.FraySolC_list[time-1])/3600.)
             self.FrayTer_list.append((self.FrayTerC_list[time]-self.FrayTerC_list[time-1])/3600.)

             if self.full :
# Qv tendencies => decumulation of cumulated flux
               self.qflux_trans_turb_list.append((self.qfluxC_trans_turb_list[time]-self.qfluxC_trans_turb_list[time-1])*24.)
               self.qflux_trans_conv_list.append((self.qfluxC_trans_conv_list[time]-self.qfluxC_trans_conv_list[time-1])*24.)
               self.qflux_cond_conv_l_list.append((self.qfluxC_cond_conv_l_list[time]-self.qfluxC_cond_conv_l_list[time-1])*24.)
               self.qflux_cond_conv_i_list.append((self.qfluxC_cond_conv_i_list[time]-self.qfluxC_cond_conv_i_list[time-1])*24.)
               self.qflux_cond_stra_l_list.append((self.qfluxC_cond_stra_l_list[time]-self.qfluxC_cond_stra_l_list[time-1])*24.)
               self.qflux_cond_stra_i_list.append((self.qfluxC_cond_stra_i_list[time]-self.qfluxC_cond_stra_i_list[time-1])*24.)
               self.qflux_evap_stra_l_list.append((self.qfluxC_evap_stra_l_list[time]-self.qfluxC_evap_stra_l_list[time-1])*24.)
               self.qflux_evap_stra_i_list.append((self.qfluxC_evap_stra_i_list[time]-self.qfluxC_evap_stra_i_list[time-1])*24.)

# H tendencies
               self.hflux_solar_list.append((self.hfluxC_solar_list[time]-self.hfluxC_solar_list[time-1])*24.)
               self.hflux_therm_list.append((self.hfluxC_therm_list[time]-self.hfluxC_therm_list[time-1])*24.)
               self.hflux_turb_list.append((self.hfluxC_turb_list[time]-self.hfluxC_turb_list[time-1])*24.)
               self.hflux_conv_list.append((self.hfluxC_conv_list[time]-self.hfluxC_conv_list[time-1])*24.)

# U and V tendencies  
               self.Uflux_conv_list.append((self.UfluxC_conv_list[time]-self.UfluxC_conv_list[time-1])*24.)
               self.Vflux_conv_list.append((self.VfluxC_conv_list[time]-self.VfluxC_conv_list[time-1])*24.)
               self.Uflux_turb_list.append((self.UfluxC_turb_list[time]-self.UfluxC_turb_list[time-1])*24.)
               self.Vflux_turb_list.append((self.VfluxC_turb_list[time]-self.VfluxC_turb_list[time-1])*24.)
               self.Uflux_gwd_list.append((self.UfluxC_gwd_list[time]-self.UfluxC_gwd_list[time-1])*24.)
               self.Vflux_gwd_list.append((self.VfluxC_gwd_list[time]-self.VfluxC_gwd_list[time-1])*24.)
          
          
             self.Pl_list.append((self.Pacc_list[time]-self.Pacc_list[time-1])*24.) 
             self.Psnow_list.append((self.Psnowacc_list[time]-self.Psnowacc_list[time-1])*24.) 
             if self.l_tctresi: self.tctresi_list.append((self.tctresiC_list[time]-self.tctresiC_list[time-1])/3600.)
             if self.l_tqvresi: self.tqvresi_list.append((self.tqvresiC_list[time]-self.tqvresiC_list[time-1])/3600.)
             if self.l_tuuresi: self.tuuresi_list.append((self.tuuresiC_list[time]-self.tuuresiC_list[time-1])/3600.)
             if self.l_tvvresi: self.tvvresi_list.append((self.tvvresiC_list[time]-self.tvvresiC_list[time-1])/3600.)
             time=time+1
          if self.l_tctresi: self.tctresi_list.append((self.tctresiC_list[time-1]-self.tctresiC_list[time-2])/3600.)
          if self.l_tqvresi: self.tqvresi_list.append((self.tqvresiC_list[time-1]-self.tqvresiC_list[time-2])/3600.)
          if self.l_tuuresi: self.tuuresi_list.append((self.tuuresiC_list[time-1]-self.tuuresiC_list[time-2])/3600.)
          if self.l_tvvresi: self.tvvresi_list.append((self.tvvresiC_list[time-1]-self.tvvresiC_list[time-2])/3600.)

           
        ilev = 1   
        while ilev <= self.nlevp1:
           self.nlevp1_list.append(ilev)
           ilev=ilev+1
        ilev = 1   
        while ilev <= self.nlev:
           self.nlev_list.append(ilev)
           ilev=ilev+1


# Calcul de la Tsrad au pas demi-horaire pour ARPEGE/AROME a partir du rayonnement

        if self.nsurf == 1:

          tsrad = []
          time = 0
          while time < self.time_slot:
             #t_soil.append(np.sqrt(np.sqrt((self.LWd_list[time] - self.LWn_list[time])/(self.epsilon*self.sigma))))
             tempo = (self.LWu_list[time] - (1. - self.epsilon)*self.LWd_list[time])/(self.epsilon*self.sigma)
             tsrad.append(np.sqrt(np.sqrt(max(0,tempo))))
             time=time+1
# Interpolation sur les slot horaires  
          self.tsrad_h.append(tsrad[0])
          time = 1     
          while time < self.time_slot:   
             self.tsrad_h.append(0.5*(tsrad[time-1]+tsrad[time]))
             time=time+1
        
          self.tsrad_h.append(tsrad[self.time_slot-1]) 

        
# Calcul des A et B
        self.compute_A_et_B(ddh_list,self.tmax)
        
        
#-----------------------------------------------------------------------------------#        
   def read_lfa_file_main(self,ddh_list,time,time_indicator):
     
        ddh_file = ddh_list[time]
        q_unit = LFA4py.wlfaouv(ddh_file,'R')

        DeltaP_sur_g = self.read_lfa_field(q_unit,'VPP'+time_indicator,1)
        self.DeltaP_sur_g_list.append(DeltaP_sur_g)
        
# Calcul de P_half et P_full        
        DeltaP = DeltaP_sur_g*self.Rg
        self.DeltaP_list.append(DeltaP)

        P_half_list = [0]
        P_full_list = []
        
        ilev=1
        while ilev < self.nlevp1:
          P_half_list.append(reduce(add, DeltaP[:ilev]))       
          P_full_list.append(P_half_list[ilev] - DeltaP[ilev-1]/2.)
          ilev=ilev+1
          
        DeltaP_sur_P =  DeltaP / P_full_list
        self.P_half_list.append(P_half_list) 
        self.P_full_list.append(P_full_list)  
                  
#  Time of slot
        echeance = self.read_lfa_field(q_unit,'ECHEANCE') 

        if time_indicator == '1': 
           self.time_list.append(int(self.read_lfa_field(q_unit,'ECHEANCE')))
           self.timef_list.append(int(self.read_lfa_field(q_unit,'ECHEANCE'))-1800)


        vuu = self.read_lfa_field(q_unit,'VUU'+time_indicator,1) 
        vvv = self.read_lfa_field(q_unit,'VVV'+time_indicator,1) 

        self.vuui_list.append(vuu)
        self.vvvi_list.append(vvv)

        vuu = vuu / DeltaP_sur_g
        vvv = vvv / DeltaP_sur_g

        ddff = self.ddff_from_u_and_v(np.array(vuu),np.array(vvv))
        self.dd_list.append(ddff[1])
        self.ff_list.append(ddff[0])
    
        self.vuu_list.append(vuu)
        self.vvv_list.append(vvv)
        
        vqvi = self.read_lfa_field(q_unit,'VQV'+time_indicator,1)
        vqli = self.read_lfa_field(q_unit,'VQL'+time_indicator,1)
        vhri = self.read_lfa_field(q_unit,'VHR'+time_indicator,1)    ###mod
        self.vqvi_list.append(vqvi)
        self.vqli_list.append(vqli)
        self.vhri_list.append(vhri)
        self.vqv_list.append(vqvi / DeltaP_sur_g)
        self.vql_list.append(vqli / DeltaP_sur_g)
        self.vhr_list.append(vhri / DeltaP_sur_g)

        if self.nsurf != 1: 
           vqni = self.read_lfa_field(q_unit,'VQI'+time_indicator,1)
        else:   
           vqni = self.read_lfa_field(q_unit,'VQN'+time_indicator,1)
           
        self.vqni_list.append(vqni)   
        self.vqn_list.append(vqni / DeltaP_sur_g)

        vnti = self.read_lfa_field(q_unit,'VNT'+time_indicator,1)
        vep  = self.read_lfa_field(q_unit,'VEP'+time_indicator,1)
        vomi = self.read_lfa_field(q_unit,'VOM'+time_indicator,1)

        try:
           #self.t_soil_h.append(self.read_lfa_field(q_unit,'SVTS'+time_indicator,1))
           self.tsurf_list.append(self.read_lfa_field(q_unit,'SVTS'+time_indicator,1))
        except RuntimeError:
           if self.debug : print('SVTS not found in file')
        except NameError:
           if self.debug : print('SVTS not found in file')
        
        self.vnti_list.append(vnti)
        self.vep_list.append(vep)
        self.vomi_list.append(vomi)
        self.vnt_list.append(vnti / DeltaP_sur_g)
        self.vzf_list.append(vep  / DeltaP)
        self.vom_list.append(vomi / DeltaP_sur_g)

        if self.nsurf == 1:
          self.Phi_surf_list.append(self.read_lfa_field(q_unit,'SVGFS05',1))
          self.lsm_list.append(self.read_lfa_field(q_unit,'S01_0',1))
          self.Snow_list.append(self.read_lfa_field(q_unit,'S06_'+time_indicator,1) / 1000.)

          if time_indicator == '1':
            self.gz0d_list.append(self.read_lfa_field(q_unit,'SVGFS06',1) / self.Rg)
            self.gz0t_list.append(self.read_lfa_field(q_unit,'SVGFS07',1) / self.Rg)
            self.pblh_list.append(self.read_lfa_field(q_unit,'SVGFS09',1))
            self.tsurf_list.append(self.read_lfa_field(q_unit,'SVGFS10',1))
            self.t2m_list.append(self.read_lfa_field(q_unit,'SVGFS01',1))
            self.q2m_list.append(self.read_lfa_field(q_unit,'SVGFS02',1))
            #self.hu2m_list.append(self.read_lfa_field(q_unit,'SVHUCLS',1))
            self.u10m_list.append(self.read_lfa_field(q_unit,'SVGFS03',1))
            self.v10m_list.append(self.read_lfa_field(q_unit,'SVGFS04',1))
            self.SWdC_list.append(self.read_lfa_field(q_unit,'SFGFS01',1))
            self.LWdC_list.append(self.read_lfa_field(q_unit,'SFGFS02',1))
            self.albedo_list.append(self.read_lfa_field(q_unit,'SVGFS08',1))
#            self.SWnC_list.append(self.read_lfa_field(q_unit,'G01',1))
            self.FraySolC_list.append(self.read_lfa_field(q_unit,'FCTRAYSOL'+time_indicator,1))
            self.SWnC_list.append(self.FraySolC_list[time][self.nlev])
#            self.LWnC_list.append(self.read_lfa_field(q_unit,'G02',1))
            self.FrayTerC_list.append(self.read_lfa_field(q_unit,'FCTRAYTER'+time_indicator,1))
            self.LWnC_list.append(self.FrayTerC_list[time][self.nlev])
            self.LheatC_list.append(self.read_lfa_field(q_unit,'SFGFS04',1)+self.read_lfa_field(q_unit,'SFGFS05',1))
            self.SheatC_list.append(self.read_lfa_field(q_unit,'SFGFS03',1))
            self.FpreciplC_list.append(self.read_lfa_field(q_unit,'FQTPRECISTL',1)+self.read_lfa_field(q_unit,'FQTPRECICOL',1))
            self.FprecipsC_list.append(self.read_lfa_field(q_unit,'FQTPRECISTN',1)+self.read_lfa_field(q_unit,'FQTPRECICON',1))
            if self.full :
              self.FpreciplconC_list.append(self.read_lfa_field(q_unit,'FQTPRECICOL',1))
              self.FprecipsconC_list.append(self.read_lfa_field(q_unit,'FQTPRECICON',1))
              self.FpreciplstrC_list.append(self.read_lfa_field(q_unit,'FQTPRECISTL',1))
              self.FprecipsstrC_list.append(self.read_lfa_field(q_unit,'FQTPRECISTN',1))
#            self.Pacc_list.append(self.read_lfa_field(q_unit,'G08',1)+self.read_lfa_field(q_unit,'G10',1))
            self.Pacc_list.append(self.FpreciplC_list[time][self.nlev])
#            self.Psnowacc_list.append(self.read_lfa_field(q_unit,'G09',1)+self.read_lfa_field(q_unit,'G11',1))
            self.Psnowacc_list.append(self.FprecipsC_list[time][self.nlev])
            self.Ustressacc_list.append(self.read_lfa_field(q_unit,'FUUTUR',1)[self.nlev])
            self.Vstressacc_list.append(self.read_lfa_field(q_unit,'FVVTUR',1)[self.nlev])

            if self.full :            
# Qv tendencies => cumulated flux, as red in ddh files

              self.qfluxC_trans_turb_list.append(self.read_lfa_field(q_unit,'FQVTUR',1))
              self.qfluxC_trans_conv_list.append(self.read_lfa_field(q_unit,'FQVTURCONV',1))
              self.qfluxC_cond_conv_l_list.append(self.read_lfa_field(q_unit,'FQTCONDECOL',1))
              self.qfluxC_cond_conv_i_list.append(self.read_lfa_field(q_unit,'FQTCONDECON',1))
              self.qfluxC_cond_stra_l_list.append(self.read_lfa_field(q_unit,'FQTCONDESTL',1))
              self.qfluxC_cond_stra_i_list.append(self.read_lfa_field(q_unit,'FQTCONDESTN',1))
              self.qfluxC_evap_stra_l_list.append(self.read_lfa_field(q_unit,'FQVEVPL',1))
              self.qfluxC_evap_stra_i_list.append(self.read_lfa_field(q_unit,'FQVEVPN',1))
            
# H tendencies
        
              self.hfluxC_solar_list.append(self.read_lfa_field(q_unit,'FCTRAYSOL1',1))
              self.hfluxC_therm_list.append(self.read_lfa_field(q_unit,'FCTRAYTER1',1))
              self.hfluxC_turb_list.append(self.read_lfa_field(q_unit,'FCTTUR',1))
              self.hfluxC_conv_list.append(self.read_lfa_field(q_unit,'FCTTURCONV',1))
          
#  U and V tendencies

              self.UfluxC_conv_list.append(self.read_lfa_field(q_unit,'FUUTURCONV',1))
              self.VfluxC_conv_list.append(self.read_lfa_field(q_unit,'FVVTURCONV',1))
              self.UfluxC_turb_list.append(self.read_lfa_field(q_unit,'FUUTUR',1))
              self.VfluxC_turb_list.append(self.read_lfa_field(q_unit,'FVVTUR',1))
              self.UfluxC_gwd_list.append(self.read_lfa_field(q_unit,'FUUONDEGREL',1))
              self.VfluxC_gwd_list.append(self.read_lfa_field(q_unit,'FVVONDEGREL',1))
            
             
            if self.l_tqvresi: self.tqvresiC_list.append(self.read_lfa_field(q_unit,'TQVRESI',1) * (echeance/86400.)/1000.)
            if self.l_tuuresi: self.tuuresiC_list.append(self.read_lfa_field(q_unit,'TUURESI',1) * echeance/86400.)
            if self.l_tvvresi: self.tvvresiC_list.append(self.read_lfa_field(q_unit,'TVVRESI',1) * echeance/86400.)
        
        elif self.nsurf == 2:    
          self.Phi_surf_list.append(self.read_lfa_field(q_unit,'SVOROG',1))
          self.lsm_list.append(self.read_lfa_field(q_unit,'SVLSM0',1))
          self.Snow_list.append(self.read_lfa_field(q_unit,'SVWN'+ time_indicator,1) / 1000.)

          if time_indicator == '1':
            self.gz0d_list.append(self.read_lfa_field(q_unit,'SVGZ01',1) / self.Rg)
            self.gz0t_list.append(self.read_lfa_field(q_unit,'SVGZH1',1) / self.Rg)
            self.pblh_list.append(self.read_lfa_field(q_unit,'SVPBLH',1))
            self.t2m_list.append(self.read_lfa_field(q_unit,'SVTCLS',1))
            self.q2m_list.append(self.read_lfa_field(q_unit,'SVQCLS',1))
            self.hu2m_list.append(self.read_lfa_field(q_unit,'SVHUCLS',1))
            self.u10m_list.append(self.read_lfa_field(q_unit,'SVUCLS',1))
            self.v10m_list.append(self.read_lfa_field(q_unit,'SVVCLS',1))
            self.SWdC_list.append(self.read_lfa_field(q_unit,'SFRAYSOLR',1))
            self.LWdC_list.append(self.read_lfa_field(q_unit,'SFRAYTHDS',1))
            self.albedo_list.append(self.read_lfa_field(q_unit,'SVALB',1))
            self.SWnC_list.append(self.read_lfa_field(q_unit,'SFRAYSODW',1))
            self.LWnC_list.append(self.read_lfa_field(q_unit,'SFRAYTHUP',1))
            self.LheatC_list.append(self.read_lfa_field(q_unit,'SFCHLATLI',1))
            self.SheatC_list.append(self.read_lfa_field(q_unit,'SFCHSENS',1))
            self.Pacc_list.append(self.read_lfa_field(q_unit,'SFPRELIGE',1))
            self.Psnowacc_list.append(self.read_lfa_field(q_unit,'SFPRENEGE',1))
            self.FpreciplC_list.append(self.read_lfa_field(q_unit,'FQTPRECISTL',1))
            self.FprecipsC_list.append(self.read_lfa_field(q_unit,'FQTPRECISTN',1))
            if self.full :
              self.FpreciplstrC_list.append(self.read_lfa_field(q_unit,'FQTPRECISTL',1))
              self.FprecipsstrC_list.append(self.read_lfa_field(q_unit,'FQTPRECISTN',1))
            self.Ustressacc_list.append(self.read_lfa_field(q_unit,'SFUUTUR',1))
            self.Vstressacc_list.append(self.read_lfa_field(q_unit,'SFVVTUR',1))
        else:
            self.Phi_surf_list.append(0.)
           

        Ps_sur_g = sum(DeltaP_sur_g)
        if time_indicator == '1':
            DeltaP_sur_g0 = self.read_lfa_field(q_unit,'VPP0',1)
            Ps_sur_g0 = sum(DeltaP_sur_g0)
            if np.fabs(Ps_sur_g0 - Ps_sur_g) > self.maxdps:
               self.maxdps = np.fabs(Ps_sur_g0 - Ps_sur_g)
               self.tmax = time

        
        self.Ps_list.append(Ps_sur_g*self.Rg)
        Cp       = self.Cpd + (self.Cpv-self.Cpd)*self.vqv_list[time]
        R        = self.Rd  + (self.Rv - self.Rd)*self.vqv_list[time]
        npR_sur_cp = np.array(R)/np.array(Cp)
        vcti = self.read_lfa_field(q_unit,'VCT'+time_indicator,1)
        self.vcti_list.append(vcti)
        self.vct_list.append(vcti / DeltaP_sur_g / Cp)
        if time_indicator == '1':
           if self.l_tctresi: self.tctresiC_list.append(self.read_lfa_field(q_unit,'TCTRESI',1) * echeance/86400.)
        npT = np.array(self.vct_list[time])
        npTheta  = npT*(self.P0/np.array(self.P_full_list[time]))**npR_sur_cp
        npThetaV = npTheta*np.array(R)/self.Rd
        self.theta_list.append(npTheta)
        self.thetav_list.append(npThetaV)

        RT       = R * self.vct_list[time] 
        gDeltaZ  = DeltaP_sur_P * RT
        Z_top    = (sum(gDeltaZ) + self.Phi_surf_list[time]) / self.Rg
        
        try :
          Zh       = [Z_top[0]] # <= Arpege formulation
        except IndexError:  
          Zh       = [Z_top]    # <= Arome case ... I don't know why !

        ilev=1
        while ilev < self.nlevp1:
           Zh.append(Zh[ilev-1]-gDeltaZ[ilev-1]/self.Rg)
           ilev = ilev+1
        
        self.vzh_list.append(Zh)
        LFA4py.wlfafer(q_unit)   

#-----------------------------------------------------------------------------------#        
   def write_file(self):
   
      print('Write NetCdf file ', self.file_name)
# Open file
      myncdf = netCDF4.Dataset(self.file_name, 'w', format='NETCDF4')

# Create global attribute       
      myncdf.title  = self.file_name
      myncdf.dataID = self.point_name
      myncdf.source = 'Initial conditions from ' + self.model + ' forecast.'
      myncdf.creator = 'lfa2nc ' + self.version + ' (Yves Bouteloup, Météo-France)'
      myncdf.details = self.model + ' ' + str(self.nlev) +  ' levels model'
      myncdf.distance = '{:.1f}'.format(self.dist)
      ddt = datetime.now()
      myncdf.NetCdf_creation_date = ddt.strftime("%A, %d. %B %Y %I:%M%p")
      myncdf.coor_par_a = np.array(self.A_list[:],dtype='f4')
      myncdf.coor_par_b = np.array(self.B_list[:],dtype='f4')
      
# Create dimension            
      time_dim   = myncdf.createDimension('time'  , self.time_slot+1) 
      time_f_dim  = myncdf.createDimension('time_f' , self.time_slot) 
      level      = myncdf.createDimension('nlev', self.nlev)
      level1     = myncdf.createDimension('nlevp1', self.nlevp1)
      levels     = myncdf.createDimension('nlevs', 1)
# Create variables   
      nlev         = myncdf.createVariable('nlev',np.int32  , ('nlev',))
      nlev.long_name = 'Model full levels' ; nlev.units   = 'count'  
      nlevp1       = myncdf.createVariable('nlevp1',np.int32  , ('nlevp1',))
      nlevp1.long_name = 'Model half levels' ; nlevp1.units   = 'count'  
      date         = myncdf.createVariable('date',np.int32  , ('time',), fill_value=-9999)
      date.long_name = 'Date'; date.units   = 'yyymmdd' ; date.associate = 'tim'    
      time         = myncdf.createVariable('time',np.int32, ('time',), fill_value=-9999)
      time.associate = 'tim';time.units   = self.time_unit ; time.long_name = 'Time'                       
      time_f        = myncdf.createVariable('time_f',np.int32, ('time_f',), fill_value=-9999)
      time_f.associate = 'tim_f';time_f.units   = self.time_unit ; time_f.long_name = 'Flux Time'
      second       = myncdf.createVariable('second',np.int32, ('time',),fill_value=-9999)
      second.long_name = 'Second';second.units = 's' 
      lat          = myncdf.createVariable('lat',np.float32, ('time',))
      lat.units    = 'deg N'   ; lat.longname = 'Latitude'
      lon          = myncdf.createVariable('lon',np.float32, ('time',))
      lon.units    = 'deg E'   ; lon.longname = 'Longitude'
      latr         = myncdf.createVariable('latr',np.float32, ('time',))
      latr.units   = 'deg N'   ; latr.longname = 'Latitude of request'
      lonr         = myncdf.createVariable('lonr',np.float32, ('time',))
      lonr.units   = 'deg E'   ; lonr.longname = 'Longitude of request'
      Z_full       = myncdf.createVariable('height_f',np.float32, ('time','nlev',),fill_value=-9999)
      Z_full.units = 'm'    ; Z_full.longname = 'Height - full levels'
      Z_half       = myncdf.createVariable('height_h',np.float32, ('time','nlevp1',),fill_value=-9999)
      Z_half.units = 'm'    ; Z_half.longname = 'Height - half levels'
      P_full       = myncdf.createVariable('pressure_f',np.float32, ('time','nlev',),fill_value=-9999)
      P_full.units = 'Pa'    ; P_full.longname = 'Pressure - full levels'
      P_half       = myncdf.createVariable('pressure_h',np.float32, ('time','nlevp1',),fill_value=-9999)
      P_half.units = 'Pa'    ; P_half.longname = 'Pressure - half levels'
      
      if self.nsurf != 0:
        h_soil       = myncdf.createVariable('h_soil',np.float32, ('nlevs',))
        h_soil.units    = 'm'   ; h_soil.longname = 'Soil layer thickness'
        t_soil       = myncdf.createVariable('t_soil',np.float32, ('time','nlevs',))
        t_soil.units    = 'K'   ; t_soil.longname = 'Soil Temperature'
        q_soil       = myncdf.createVariable('q_soil',np.float32, ('time','nlevs',))
        q_soil.units    = 'm^3/m^3'   ; q_soil.longname = 'Soil Moisture'
        snow         = myncdf.createVariable('snow',np.float32, ('time',))
        snow.units   = 'm, liquid equivalent'   ; snow.longname = 'Snow Depth'

      Ps           = myncdf.createVariable('ps',np.float32, ('time',))
      Ps.units     = 'Pascal'  ; Ps.longname = 'Surface Pressure'
      U            = myncdf.createVariable('u',np.float32, ('time','nlev',))
      U.units      = 'm/s'    ; U.longname = 'Zonal Wind'
      V            = myncdf.createVariable('v',np.float32, ('time','nlev',))
      V.units      = 'm/s'    ; V.longname = 'Meridional Wind'
      Vdir         = myncdf.createVariable('Vdir',np.float32, ('time','nlev',))
      Vdir.units   = 'degrees from N'    ; Vdir.longname = 'Wind direction'
      Vamp         = myncdf.createVariable('Vamp',np.float32, ('time','nlev',))
      Vamp.units   = 'm/s'    ; Vamp.longname = 'Wind speed'
      UG           = myncdf.createVariable('ug',np.float32, ('time','nlev',))
      UG.units     = 'm/s'    ; UG.longname = 'Geostrophic U Wind'
      VG           = myncdf.createVariable('vg',np.float32, ('time','nlev',))
      VG.units     = 'm/s'    ; VG.longname = 'Geostrophic V Wind'
      T            = myncdf.createVariable('t',np.float32, ('time','nlev',))
      T.units      = 'K'       ; T.longname = 'Temperature'
      Theta        = myncdf.createVariable('theta',np.float32, ('time','nlev',))
      Theta.units  = 'K'       ; Theta.longname = 'Potential Temperature'
      Thv          = myncdf.createVariable('thv',np.float32, ('time','nlev',))
      Thv.units    = 'K'       ; Thv.longname = 'Virtual Potential Temperature'
      Q            = myncdf.createVariable('q',np.float32, ('time','nlev',))
      Q.units      = 'kg/kg'  ; Q.longname = 'Water Vapor Mixing Ratio'
      Qv           = myncdf.createVariable('qv',np.float32, ('time','nlev',))
      Qv.units     = 'kg/kg'  ; Qv.longname = 'Water Vapor Mixing Ratio'
      Ql           = myncdf.createVariable('ql',np.float32, ('time','nlev',))
      Ql.units     = 'kg/kg'  ; Ql.longname = 'Liquid Water Mixing Ratio'
      Hr           = myncdf.createVariable('hr',np.float32, ('time','nlev',))
      Hr.units     = '0-1'  ; Hr.longname = 'relativ humidity'
      Qi           = myncdf.createVariable('qi',np.float32, ('time','nlev',))
      Qi.units     = 'kg/kg'  ; Qi.longname = 'Ice Water Mixing Ratio'
      Cf           = myncdf.createVariable('cloud_fraction',np.float32, ('time','nlev',))
      Cf.units     = '0-100'     ; Cf.longname = 'Cloud Fraction'
      Cc           = myncdf.createVariable('cc',np.float32, ('time','nlev',))
      Cc.units     = '0-1'     ; Cc.longname = 'Cloud Fraction'
      omega        = myncdf.createVariable('omega',np.float32, ('time','nlev',))
      omega.units  = 'Pa/s'  ; omega.longname = 'Vertical Pressure Velocity'


      if self.full :     
# Qv tendencies 
             
        qtend_trans_turb = myncdf.createVariable('qtend_trans_turb',np.float32, ('time_f','nlev',))
        qtend_trans_turb.units     = 'kg/kg/day'     ; qtend_trans_turb.longname = 'Turbulent transport of Qv'
        qtend_trans_conv = myncdf.createVariable('qtend_trans_conv',np.float32, ('time_f','nlev',))
        qtend_trans_conv.units     = 'kg/kg/day'     ; qtend_trans_conv.longname = 'Convective transport of Qv'
        qtend_cond_conv_l = myncdf.createVariable('qtend_cond_conv_l',np.float32, ('time_f','nlev',))
        qtend_cond_conv_l.units     = 'kg/kg/day'     ; qtend_cond_conv_l.longname = 'Convective condensation liquid'
        qtend_cond_conv_i = myncdf.createVariable('qtend_cond_conv_i',np.float32, ('time_f','nlev',))
        qtend_cond_conv_i.units     = 'kg/kg/day'     ; qtend_cond_conv_i.longname = 'Convective condensation solid'
        qtend_cond_conv_t = myncdf.createVariable('qtend_cond_conv_t',np.float32, ('time_f','nlev',))
        qtend_cond_conv_t.units     = 'kg/kg/day'     ; qtend_cond_conv_t.longname = 'Total convective condensation'
        qtend_cond_stra_l = myncdf.createVariable('qtend_cond_stra_l',np.float32, ('time_f','nlev',))
        qtend_cond_stra_l.units     = 'kg/kg/day'     ; qtend_cond_stra_l.longname = 'Stratiforme condensation liquid'
        qtend_cond_stra_i = myncdf.createVariable('qtend_cond_stra_i',np.float32, ('time_f','nlev',))
        qtend_cond_stra_i.units     = 'kg/kg/day'     ; qtend_cond_stra_i.longname = 'Stratiforme condensation solid'
        qtend_cond_stra_t = myncdf.createVariable('qtend_cond_stra_t',np.float32, ('time_f','nlev',))
        qtend_cond_stra_t.units     = 'kg/kg/day'     ; qtend_cond_stra_t.longname = 'Stratiforme total condensation'
        qtend_evap_stra_l = myncdf.createVariable('qtend_evap_stra_l',np.float32, ('time_f','nlev',))
        qtend_evap_stra_l.units     = 'kg/kg/day'     ; qtend_evap_stra_l.longname = 'Stratiforme evaporation liquid'
        qtend_evap_stra_i = myncdf.createVariable('qtend_evap_stra_i',np.float32, ('time_f','nlev',))
        qtend_evap_stra_i.units     = 'kg/kg/day'     ; qtend_evap_stra_i.longname = 'Stratiforme evaporation solid'
        qtend_evap_stra_t = myncdf.createVariable('qtend_evap_stra_t',np.float32, ('time_f','nlev',))
        qtend_evap_stra_t.units     = 'kg/kg/day'     ; qtend_evap_stra_t.longname = 'Stratiforme total evaporation'
        qtend_phy_total = myncdf.createVariable('qtend_phy_total',np.float32, ('time_f','nlev',))
        qtend_phy_total.units     = 'kg/kg/day'     ; qtend_phy_total.longname = 'Qv total physical tendencies'
        qtend_total = myncdf.createVariable('qtend_total',np.float32, ('time_f','nlev',))
        qtend_total.units     = 'kg/kg/day'     ; qtend_total.longname = 'Qv total tendencies'
        qtend_dyn = myncdf.createVariable('qtend_dyn',np.float32, ('time_f','nlev',))
        qtend_dyn.units     = 'kg/kg/day'     ; qtend_dyn.longname = 'Qv dynamical tendencies'

# H tendencies

        htend_solar = myncdf.createVariable('htend_solar',np.float32, ('time_f','nlev',))
        htend_solar.units     = 'J/kg/day'     ; htend_solar.longname = 'Enthalpy solar tendencies'
        htend_therm = myncdf.createVariable('htend_therm',np.float32, ('time_f','nlev',))
        htend_therm.units     = 'J/kg/day'     ; htend_therm.longname = 'Enthalpy thermal tendencies'
        htend_turb = myncdf.createVariable('htend_turb',np.float32, ('time_f','nlev',))
        htend_turb.units     = 'J/kg/day'     ; htend_turb.longname = 'Enthalpy turbulent tendencies'
        htend_conv = myncdf.createVariable('htend_conv',np.float32, ('time_f','nlev',))
        htend_conv.units     = 'J/kg/day'     ; htend_conv.longname = 'Enthalpy convective tendencies'
        htend_cond_l = myncdf.createVariable('htend_cond_l',np.float32, ('time_f','nlev',))
        htend_cond_l.units     = 'J/kg/day'     ; htend_cond_l.longname = 'Liquid-vapor phase change Enthalpy tendencies'
        htend_cond_i = myncdf.createVariable('htend_cond_i',np.float32, ('time_f','nlev',))
        htend_cond_i.units     = 'J/kg/day'     ; htend_cond_i.longname = 'Ice-vapor phase change Enthalpy tendencies'
        htend_rain = myncdf.createVariable('htend_rain',np.float32, ('time_f','nlev',))
        htend_rain.units     = 'J/kg/day'     ; htend_rain.longname = 'Enthalpy tendencies due to falling of rain'
        htend_snow = myncdf.createVariable('htend_snow',np.float32, ('time_f','nlev',))
        htend_snow.units     = 'J/kg/day'     ; htend_snow.longname = 'Enthalpy tendencies due to falling of snow'
        htend_total = myncdf.createVariable('htend_total',np.float32, ('time_f','nlev',))
        htend_total.units     = 'J/kg/day'     ; htend_total.longname = 'Total Enthalpy tendencies'
        htend_phy_total = myncdf.createVariable('htend_phy_total',np.float32, ('time_f','nlev',))
        htend_phy_total.units     = 'J/kg/day'     ; htend_phy_total.longname = 'Enthalpy total physical tendencies'
        htend_dyn = myncdf.createVariable('htend_dyn',np.float32, ('time_f','nlev',))
        htend_dyn.units     = 'J/kg/day'     ; htend_dyn.longname = 'Enthalpy dynamical tendencies'
        htend_dyn_ray = myncdf.createVariable('htend_dyn_ray',np.float32, ('time_f','nlev',))
        htend_dyn_ray.units     = 'J/kg/day'     ; htend_dyn_ray.longname = 'Enthalpy dynamical and radiative tendencies'

# T tendencies

        ttend_dyn = myncdf.createVariable('ttend_dyn',np.float32, ('time_f','nlev',))
        ttend_dyn.units     = 'K/day'     ; ttend_dyn.longname = 'Temperature dynamical tendencies'
        ttend_dyn_ray = myncdf.createVariable('ttend_dyn_ray',np.float32, ('time_f','nlev',))
        ttend_dyn_ray.units     = 'K/day'     ; ttend_dyn_ray.longname = 'Temperature dynamical and radiative tendencies'

# U and V tendencies  

        utend_conv = myncdf.createVariable('utend_conv',np.float32, ('time_f','nlev',))
        utend_conv.units     = 'm/s/day'     ; utend_conv.longname = 'Zonal wind convective tendencies'
        vtend_conv = myncdf.createVariable('vtend_conv',np.float32, ('time_f','nlev',))
        vtend_conv.units     = 'm/s/day'     ; vtend_conv.longname = 'Meridional wind convective tendencies'
        utend_turb = myncdf.createVariable('utend_turb',np.float32, ('time_f','nlev',))
        utend_turb.units     = 'm/s/day'     ; utend_turb.longname = 'Zonal wind turbulent tendencies'
        vtend_turb = myncdf.createVariable('vtend_turb',np.float32, ('time_f','nlev',))
        vtend_turb.units     = 'm/s/day'     ; vtend_turb.longname = 'Meridional wind turbulent tendencies'
        utend_gwd = myncdf.createVariable('utend_gwd',np.float32, ('time_f','nlev',))
        utend_gwd.units     = 'm/s/day'     ; utend_gwd.longname = 'Zonal wind GWD tendencies'
        vtend_gwd = myncdf.createVariable('vtend_gwd',np.float32, ('time_f','nlev',))
        vtend_gwd.units     = 'm/s/day'     ; vtend_gwd.longname = 'Meridional wind GWD tendencies'
        utend_phy_total = myncdf.createVariable('utend_phy_total',np.float32, ('time_f','nlev',))
        utend_phy_total.units     = 'm/s/day'     ; utend_phy_total.longname = 'Zonal wind total physical tendencies'
        vtend_phy_total = myncdf.createVariable('vtend_phy_total',np.float32, ('time_f','nlev',))
        vtend_phy_total.units     = 'm/s/day'     ; vtend_phy_total.longname = 'Merifional wind total physical tendencies'
        utend_total = myncdf.createVariable('utend_total',np.float32, ('time_f','nlev',))
        utend_total.units     = 'm/s/day'     ; utend_total.longname = 'Zonal wind total tendencies'
        vtend_total = myncdf.createVariable('vtend_total',np.float32, ('time_f','nlev',))
        vtend_total.units     = 'm/s/day'     ; vtend_total.longname = 'Merifional wind total tendencies'
        utend_dyn = myncdf.createVariable('utend_dyn',np.float32, ('time_f','nlev',))
        utend_dyn.units     = 'm/s/day'     ; utend_dyn.longname = 'Zonal wind dynamical tendencies'
        vtend_dyn = myncdf.createVariable('vtend_dyn',np.float32, ('time_f','nlev',))
        vtend_dyn.units     = 'm/s/day'     ; vtend_dyn.longname = 'Merifional wind dynamical tendencies'

      if self.nsurf != 0:
        P            = myncdf.createVariable('P',np.float32, ('time_f'))
        P.units      = 'mm/day'  ; P.longname = 'Surface Precipitation Rate'
        Pacc         = myncdf.createVariable('Pacc',np.float32, ('time_f'))
        Pacc.units   = 'mm'      ; Pacc.longname = 'Accumulated P'
        Psnow        = myncdf.createVariable('Psnow',np.float32, ('time_f'))
        Psnow.units  = 'mm/day'  ; Psnow.longname = 'Surface Snowfall Rate'
        Psnowacc     = myncdf.createVariable('Psnowacc',np.float32, ('time_f'))
        Psnowacc.units  = 'mm'  ; Psnowacc.longname = 'Accumulated Psnow'
        Zgs          = myncdf.createVariable('orog',np.int32, ('time',))
        Zgs.units    = 'm2/s2'   ; Zgs.longname = 'Orography (surface geopotential)'
        albedo       = myncdf.createVariable('albedo',np.float32, ('time_f',))
        albedo.units = '1'     ; albedo.longname = 'Albedo'
        mom_rough    = myncdf.createVariable('mom_rough',np.float32, ('time_f',))
        mom_rough.units = 'm'     ; mom_rough.longname = 'Momentum Roughness Length'
        heat_rough    = myncdf.createVariable('heat_rough',np.float32, ('time_f',))
        heat_rough.units = 'm'     ; heat_rough.longname = 'Heat Roughness Length'
        lsm          = myncdf.createVariable('lsm',np.float32, ('time',))
        lsm.units      = '0-1'     ; lsm.longname = 'Land Sea mask'
      
     #  if self.nsurf == 1:
        Pflux        = myncdf.createVariable('Pflux',np.float32, ('time_f','nlevp1',))
        Pflux.units  = 'mm/day'  ; Pflux.longname = 'Total Precipitation Flux'
        Pfluxs       = myncdf.createVariable('Pfluxs',np.float32, ('time_f','nlevp1',))
        Pfluxs.units = 'mm/day'  ; Pfluxs.longname = 'Total Snow Flux'

        if self.full :
          Pfluxcon     = myncdf.createVariable('Pfluxcon',np.float32, ('time_f','nlevp1',))
          Pfluxcon.units  = 'mm/day'  ; Pfluxcon.longname = 'Convective Precipitation Flux'
          Pfluxscon       = myncdf.createVariable('Pfluxscon',np.float32, ('time_f','nlevp1',))
          Pfluxscon.units = 'mm/day'  ; Pfluxscon.longname = 'Convective Snow Flux'

          Pfluxstr     = myncdf.createVariable('Pfluxstr',np.float32, ('time_f','nlevp1',))
          Pfluxstr.units  = 'mm/day'  ; Pfluxstr.longname = 'Stratiform Precipitation Flux'
          Pfluxsstr       = myncdf.createVariable('Pfluxsstr',np.float32, ('time_f','nlevp1',))
          Pfluxsstr.units = 'mm/day'  ; Pfluxsstr.longname = 'Stratiform Snow Flux'

        ustress      = myncdf.createVariable('ustress',np.float32, ('time_f'))
        ustress.units = 'm2/s2'   ; ustress.longname = 'Surface U stress'
        vstress      = myncdf.createVariable('vstress',np.float32, ('time_f'))
        vstress.units = 'm2/s2'   ; vstress.longname = 'Surface V stress'
      
      uadv      = myncdf.createVariable('uadv',np.float32, ('time','nlev',))
      if self.l_tuuresi: 
         uadv.units = 'm/s/s'
         uadv.longname = 'Zonal wind total advection (hor + vert)'
      else:
         uadv.units = ' '
         uadv.longname = 'Always 0, just for testbed'
      vadv      = myncdf.createVariable('vadv',np.float32, ('time','nlev',))
      if self.l_tvvresi: 
         vadv.units = 'm/s/s'
         vadv.longname = 'Meridional wind total advection (hor + vert)'
      else:
         vadv.units = ' '
         vadv.longname = 'Always 0, just for testbed'
      tadv      = myncdf.createVariable('tadv',np.float32, ('time','nlev',))
      if self.l_tctresi: 
         tadv.units = 'K/s'   ; 
         tadv.longname = 'Temperature total advection (hor + vert)'
      else:   
         tadv.units = ' '   ; 
         tadv.longname = 'Always 0, just for testbed'
      qadv      = myncdf.createVariable('qadv',np.float32, ('time','nlev',))
      if self.l_tqvresi:  
         qadv.units = 'kg/kg/s'
         qadv.longname = 'Water Vapor total advection (hor + vert)'
      else:   
         qadv.units = ' '
         qadv.longname = 'Always 0, just for testbed'
         
      ladv      = myncdf.createVariable('ladv',np.float32, ('time','nlev',))
      ladv.units = ' '   ; ladv.longname = 'Always 0, just for testbed'
      iadv      = myncdf.createVariable('iadv',np.float32, ('time','nlev',))
      iadv.units = ' '   ; iadv.longname = 'Always 0, just for testbed'
      
      if self.nsurf != 0:
      # if self.nsurf == 1:
        mfsfc        = myncdf.createVariable('mfsfc',np.float32, ('time_f'))
        mfsfc.units  = 'm2/s2'   ; mfsfc.longname = 'Surface Momentum flux'
      
        pblh         = myncdf.createVariable('pblh',np.float32, ('time_f',))
        pblh.units   = 'm'       ; pblh.longname = 'Planetary Boundary Layer Height'
        tsrad        = myncdf.createVariable('tsrad',np.float32, ('time',))
        tsrad.units   = 'K'       ; tsrad.longname = 'Mean Radiative Surface Temperature'
        tsurf        = myncdf.createVariable('tsurf',np.float32, ('time',))
        tsurf.units   = 'K'       ; tsurf.longname = 'Mean Surface Temperature'
        t2m          = myncdf.createVariable('t2m',np.float32, ('time',))
        t2m.units    = 'K'       ; t2m.longname = 'Temperature at 2m'
        q2m          = myncdf.createVariable('q2m',np.float32, ('time',))
        q2m.units    = 'kg/kg' ; q2m.longname = 'Specific humidity at 2m'
        if self.nsurf == 2:
          hu2m          = myncdf.createVariable('hu2m',np.float32, ('time',))
          hu2m.units    = '0-1' ; hu2m.longname = 'Relative humidity at 2m'
        u10m         = myncdf.createVariable('u10m',np.float32, ('time',))
        u10m.units   = 'm/s'   ; u10m.longname = 'Zonal wind speed at 10m'
        v10m         = myncdf.createVariable('v10m',np.float32, ('time',))
        v10m.units   = 'm/s'   ; v10m.longname = 'Meridional wind speed at 10m'
        Vamp10m      = myncdf.createVariable('Vamp10m',np.float32, ('time',))
        Vamp10m .units = 'm/s'   ; Vamp10m .longname = 'Wind speed at 10m'
        Vdir10m      = myncdf.createVariable('Vdir10m',np.float32, ('time',))
        Vdir10m .units = 'degrees from N'   ; Vdir10m .longname = 'Wind direction at 10m'
        SWd          = myncdf.createVariable('SWd',np.float32, ('time_f',))
        SWd.units    = 'W/m2'   ; SWd.longname = 'Surface downward SW rad. flux (+ downw.)'
        LWd          = myncdf.createVariable('LWd',np.float32, ('time_f',))
        LWd.units    = 'W/m2'   ; LWd.longname = 'Surface downward LW rad. flux (+ downw.)'
        SWn          = myncdf.createVariable('SWn',np.float32, ('time_f',))
        SWn.units    = 'W/m2'   ; SWn.longname = 'Surface net SW rad. flux (+ downw.)'
        LWn          = myncdf.createVariable('LWn',np.float32, ('time_f',))
        LWn.units    = 'W/m2'   ; LWn.longname = 'Surface net LW rad. flux (+ downw.)'
        LWu          = myncdf.createVariable('LWu',np.float32, ('time_f',))
        LWu.units    = 'W/m2'   ; LWu.longname = 'Surface upward LW rad. flux'
        SWu          = myncdf.createVariable('SWu',np.float32, ('time_f',))
        SWu.units    = 'W/m2'   ; SWu.longname = 'Surface upward SW rad. flux'
        sfc_sens_flx = myncdf.createVariable('sfc_sens_flx',np.float32, ('time_f',))
        sfc_sens_flx.units = 'W/m2'   ; sfc_sens_flx.longname = 'Surface Sensible Heat Flux'
        sfc_lat_flx = myncdf.createVariable('sfc_lat_flx',np.float32, ('time_f',))
        sfc_lat_flx.units = 'W/m2'   ; sfc_lat_flx.longname = 'Surface Latent Heat Flux'

# Put data in variables
      date[:] = self.ddt
      nlev[:] = self.nlev_list[:]
      nlevp1[:] = self.nlevp1_list[:]
      time[:] = self.time_list[:]
      time_f[:] = self.timef_list[:]
      second[:] = time[:] + self.hour*3600
      lat[:]    = self.lat
      lon[:]    = self.lon
      latr[:]   = self.latr
      lonr[:]   = self.lonr
      Ps[:]     = self.Ps_list[:]
      U[:,:]    = self.vuu_list[:][:]
      V[:,:]    = self.vvv_list[:][:]
      Vdir[:,:] = self.dd_list[:][:]
      Vamp[:,:] = self.ff_list[:][:]
      T[:,:]    = self.vct_list[:][:]
      Theta[:,:]    = self.theta_list[:][:]
      Thv[:,:]    = self.thetav_list[:][:]
      Q[:,:]    = self.vqv_list[:][:]
      Qv[:,:]   = self.vqv_list[:][:]
      Ql[:,:]   = self.vql_list[:][:]
      Hr[:,:]   = self.vhr_list[:][:]

      Qi[:,:]   = self.vqn_list[:][:]
      Cf[:,:]   = np.array(self.vnt_list[:][:])*100.
      Cc[:,:]   = self.vnt_list[:][:]
      omega[:,:]= self.vom_list[:][:]
      P_full[:] = self.P_full_list[:][:]
      P_half[:] = self.P_half_list[:][:]
      Z_full[:] = self.vzf_list[:][:]
      
  #    print(Z_full)
  #    print("----------")
  #    print(len(self.vzf_list[0]))  
  #    print("----------")
  #    print(Z_half)
  #    print("----------")

   #   print(len(self.vzh_list[0]))   
      Z_half[:] = self.vzh_list[:][:]      
      
      if self.nsurf != 0:
        Zgs[:]    = self.Phi_surf_list[:]
        albedo[:] = self.albedo_list[:]
        mom_rough[:] = self.gz0d_list[:]
        lsm[:] = self.lsm_list[:]
        heat_rough[:] = self.gz0t_list[:]
        pblh[:]   = self.pblh_list[:]
        tsurf[:]   = self.tsurf_list[:]
        t2m[:]   = self.t2m_list[:]
        q2m[:]   = self.q2m_list[:]
        if self.nsurf == 2:
          hu2m[:]   = self.hu2m_list[:]
        u10m[:]   = self.u10m_list[:]
        v10m[:]   = self.v10m_list[:]
        npu10m = np.array(self.u10m_list[:])
        npv10m = np.array(self.v10m_list[:])
        ddff = self.ddff_from_u_and_v(npu10m,npv10m)
        Vamp10m[:]   = ddff[0]
        Vdir10m[:]   = ddff[1]
        SWd[:]   = self.SWd_list[:]
        LWd[:]   = self.LWd_list[:]
        SWn[:]   = self.SWn_list[:]
        LWn[:]   = self.LWn_list[:]
        LWu[:]   = self.LWu_list[:]
        SWu[:]   = self.SWu_list[:]
      
#      t_soil[:,:] = np.sqrt(np.sqrt((0.965*np.array(self.LWd_list[:]) - np.array(self.LWn_list[:]))/self.sigma))
        #t_soil[:,:] = self.t_soil_h[:]
        t_soil[:,:] = self.tsurf_list[:]
        tsrad[:] = self.tsrad_h[:]
        h_soil[:] = 0.01
        q_soil[:,:] = -999.
        snow[:] = self.Snow_list[:]
        sfc_sens_flx[:] = self.Sheat_list[:]
        sfc_lat_flx[:]  = self.Lheat_list[:]
        P[:]  = self.Pl_list[:]
        Psnow[:]  = self.Psnow_list[:]
        Pacc[:]  = self.Pacc_list[:]
        Psnowacc[:]  = self.Psnowacc_list[:]

     #  if self.nsurf == 1:
        ustress[:] = self.Ustress_list[:]
        vstress[:] = self.Vstress_list[:]
        Pflux[:,:]  = self.Fprecipl_list[:][:]
        Pfluxs[:,:] = self.Fprecips_list[:][:]

        if self.full :
          Pfluxcon[:,:]  = self.Fpreciplcon_list[:][:]
          Pfluxscon[:,:] = self.Fprecipscon_list[:][:]
          Pfluxstr[:,:]  = self.Fpreciplstr_list[:][:]
          Pfluxsstr[:,:] = self.Fprecipsstr_list[:][:]

#  Qv tendencies ...        

          qtend_trans_turb[:,:]  = self.flux_2_tend(self.qflux_trans_turb_list)
          qtend_trans_conv[:,:]  = self.flux_2_tend(self.qflux_trans_conv_list)
          qtend_cond_conv_l[:,:] = self.flux_2_tend(self.qflux_cond_conv_l_list)
          qtend_cond_conv_i[:,:] = self.flux_2_tend(self.qflux_cond_conv_i_list)
          qtend_cond_conv_t[:,:] = qtend_cond_conv_l[:,:] + qtend_cond_conv_i[:,:]
          qtend_cond_stra_l[:,:] = self.flux_2_tend(self.qflux_cond_stra_l_list)
          qtend_cond_stra_i[:,:] = self.flux_2_tend(self.qflux_cond_stra_i_list)
          qtend_cond_stra_t[:,:] = qtend_cond_stra_l[:,:] + qtend_cond_stra_i[:,:]
          qtend_evap_stra_l[:,:] = self.flux_2_tend(self.qflux_evap_stra_l_list)
          qtend_evap_stra_i[:,:] = self.flux_2_tend(self.qflux_evap_stra_i_list)
          qtend_evap_stra_t[:,:] = qtend_evap_stra_l[:,:] + qtend_evap_stra_i[:,:]
          qtend_phy_total[:,:]   = qtend_trans_turb[:,:] + qtend_trans_conv[:,:] + qtend_cond_stra_t[:,:] + qtend_cond_conv_t[:,:] \
          - qtend_evap_stra_t[:,:]
          qtend_total[:,:]       = self.var_2_tend(Q)
          qtend_dyn[:,:]         = qtend_total[:,:] - qtend_phy_total[:,:]

# H tendencies

          Cp       = self.Cpd + (self.Cpv-self.Cpd)*np.array(self.vqv_list[:][:])

          htend_solar[:,:]  = self.flux_2_tend(self.hflux_solar_list)
          htend_therm[:,:]  = self.flux_2_tend(self.hflux_therm_list)
          htend_turb[:,:]   = self.flux_2_tend(self.hflux_turb_list)
          htend_conv[:,:]   = self.flux_2_tend(self.hflux_conv_list)
          htend_cond_l[:,:] = self.Rlvzer*(qtend_cond_conv_l[:,:] + qtend_cond_stra_l[:,:] - qtend_evap_stra_l[:,:])
          htend_cond_i[:,:] = self.Rlszer*(qtend_cond_conv_i[:,:] + qtend_cond_stra_i[:,:] - qtend_evap_stra_i[:,:])
         
          htend_rain[:,:]   = self.flux_2_tend((self.Rcw-self.Cpd)*Pflux[:,:]*self.var_2_varf(self.vct_list[:][:]))
          htend_snow[:,:]   = self.flux_2_tend((self.Rcs-self.Cpd)*Pfluxs[:,:]*self.var_2_varf(self.vct_list[:][:]))
        
          htend_total[:,:]       = self.var_2_tend(T)*Cp[0:self.time_slot][:]
          htend_phy_total[:,:]   = htend_solar[:,:] + htend_therm[:,:] + htend_turb[:,:] + htend_conv[:,:] - htend_cond_l[:,:] \
          - htend_cond_i[:,:] + htend_rain[:,:] + htend_snow[:,:]
          htend_dyn[:,:]         = htend_total[:,:] - htend_phy_total[:,:]
          htend_dyn_ray[:,:]     = htend_dyn[:,:] + htend_solar[:,:] + htend_therm[:,:]


# T tendencies

          ttend_dyn[:,:]     = htend_dyn[:,:] / Cp[0:self.time_slot][:]
    #    ttend_dyn[:,:]     = htend_total[:,:] / Cp[0:self.time_slot][:]
    #    ttend_dyn[:,:]     = self.var_2_tend(T)
          ttend_dyn_ray[:,:] = htend_dyn_ray[:,:] / Cp[0:self.time_slot][:]
        
# U and V tendencies

          utend_conv[:,:]  = self.flux_2_tend(self.Uflux_conv_list)
          vtend_conv[:,:]  = self.flux_2_tend(self.Vflux_conv_list)
          utend_turb[:,:]  = self.flux_2_tend(self.Uflux_turb_list)
          vtend_turb[:,:]  = self.flux_2_tend(self.Vflux_turb_list)
          utend_gwd[:,:]   = self.flux_2_tend(self.Uflux_gwd_list)
          vtend_gwd[:,:]   = self.flux_2_tend(self.Vflux_gwd_list)
        
          utend_phy_total[:,:] = utend_conv[:,:] + utend_turb[:,:] + utend_gwd[:,:]
          vtend_phy_total[:,:] = vtend_conv[:,:] + vtend_turb[:,:] + vtend_gwd[:,:]
          utend_total[:,:]     = self.var_2_tend(U)
          vtend_total[:,:]     = self.var_2_tend(V)
          utend_dyn[:,:]       = utend_total[:,:] - utend_phy_total[:,:]
          vtend_dyn[:,:]       = vtend_total[:,:] - vtend_phy_total[:,:]

# calcul de u*      
        npustress = np.array(self.Ustress_list[:])
        npvstress = np.array(self.Vstress_list[:])
        ustar     = self.ddff_from_u_and_v(npustress,npvstress)
        mfsfc[:]  = ustar[0]

# Advection, not always 0, depends on input files
      
      if self.l_tctresi: 
         tadv[:,:] = self.tctresi_list[:][:]
      else:  
         tadv[:,:] = 0. 
      if self.l_tqvresi: 
         qadv[:,:] = self.tqvresi_list[:][:]
      else:
         qadv[:,:] = 0.
      if self.l_tuuresi: 
         uadv[:,:] = self.tuuresi_list[:][:]   
      else:      
         uadv[:,:] = 0.
      if self.l_tvvresi: 
         vadv[:,:] = self.tvvresi_list[:][:]
      else:        
         vadv[:,:] = 0. 
      ladv[:,:] = 0.
      iadv[:,:] = 0.
# Geostrophic Wind (missing value -999.)   
      UG[:,:] = -999.   
      VG[:,:] = -999. 
# Close file      
      myncdf.close()

#-----------------------------------------------------------------------------------#        
   def compute_A_et_B(self,ddh_list,time):
#  def compute_A_et_B(self,ddh_list,point,time):

      ddh_file = ddh_list[time]
      q_unit = LFA4py.wlfaouv(ddh_file,'R')
   
      DeltaP0 = self.read_lfa_field(q_unit,'VPP0',1) * self.Rg   
      DeltaP1 = self.read_lfa_field(q_unit,'VPP1',1) * self.Rg   
      Ps0 = sum(DeltaP0)
      Ps1 = sum(DeltaP1)
      self.A_list  = [0.]
      self.B_list  = [0.]
      ilev = 1
      while ilev < self.nlevp1:
        P0 = reduce(add, DeltaP0[:ilev])        
        P1 = reduce(add, DeltaP1[:ilev])
        self.B_list.append(np.fabs((P1-P0)/(Ps1-Ps0)))
        self.A_list.append(max(0.,P1-self.B_list[ilev]*Ps1))
        ilev=ilev+1
        
      LFA4py.wlfafer(q_unit)   
        
#-----------------------------------------------------------------------------------#        
   def point_name_maker(self,lat,lon):
      
      ilat = int(lat*10.)
      ilon = int(lon*10.)
      if ilat==673 and ilon==266:
        place_name = 'Sodankyla'
      elif ilat==519 and ilon==49:
        place_name = 'Cabauw'
      elif ilat==487 and ilon==22:
        place_name = 'Sirta'
      elif ilat==484 and ilon==22:
        place_name = 'Sirta'
      elif ilat==435 and ilon==13:
        place_name = 'Toulouse'
      elif ilat==448 and ilon==-6:
        place_name = 'Bordeaux'
      elif ilat==448 and ilon==3594:
        place_name = 'Bordeaux'
      elif ilat==438 and ilon==44:
        place_name = 'Nimes'
      elif ilat==431 and ilon==3:
        place_name = 'Lannemezan'
      elif ilat==418 and ilon==-57:
        place_name = 'Valladolid'
      elif ilat==418 and ilon==3543:
        place_name = 'Valladolid'
      elif ilat==601 and ilon==246:
        place_name = 'Kivenlahti'
      elif ilat==626 and ilon==274:
        place_name = 'Kuopio'
      elif ilat==665 and ilon==255:
        place_name = 'Rovaniemi'
      elif ilat==-750 and ilon==1233:
        place_name = 'Dome-C'
      elif ilat==366 and ilon==-975:
        place_name = 'Arm-SGP'
      elif ilat==366 and ilon==2624:
        place_name = 'Arm-SGP'
      elif ilat==713 and ilon==-1566:
        place_name = 'Arm-Barrow'
      elif ilat==713 and ilon==2034:
        place_name = 'Arm-Barrow'
      elif ilat==-5 and ilon==1669:
        place_name = 'Arm-Nauru'
      elif ilat==-20 and ilon==1474:
        place_name = 'Arm-Manus'
      elif ilat==-124 and ilon==1308:
        place_name = 'Arm-Darwin'
      elif ilat==511 and ilon==-14:
        place_name = 'Chilbolton'
      elif ilat==511 and ilon==3586:
        place_name = 'Chilbolton'
      elif ilat==521 and ilon==141:
        place_name = 'Lindenberg'
      elif ilat==405 and ilon==157:
        place_name = 'Brienza'
      elif ilat==134 and ilon==21:
        place_name = 'Niamey'
      elif ilat==799 and ilon==-859:
        place_name = 'Eureka'
      elif ilat==799 and ilon==2741:
        place_name = 'Eureka'
      elif ilat==347 and ilon==1272:
        place_name = 'Bosung-Tower'
      elif ilat==-666 and ilon==1400:
        place_name = 'Dumont-d-Urville'
      elif ilat==-755 and ilon==-266:
        place_name = 'Halley'
      else:  
        place_name ='lat' + str(int(lat*10.)) + 'lon' + str(int(lon*10.))
        
      return place_name    
      
#-----------------------------------------------------------------------------------#        
   def ddff_from_u_and_v(self,u,v):
   
      ff = np.sqrt(u**2+v**2)
      ffcal = np.where(ff < 1.E-15, 1.E-15, ff)

      u_norm = -u/ffcal
      v_norm = -v/ffcal
      dd1 = np.arccos(v_norm)
      dd2 = 2.*np.pi - dd1     
      dd = np.where(u_norm >= 0., dd1, dd2)*180./np.pi

      return(ff,dd)

#-----------------------------------------------------------------------------------#        
   def read_lfa_field(self,q_unit,field,extract=0):
      try:
        (fieldtype, fieldlength)  = LFA4py.wlfacas(q_unit, field)
        if fieldtype == 'I4':
          (data, fieldlength) = LFA4py.wlfaleci(q_unit, field, fieldlength)  
        elif fieldtype == 'R4':   
          (data, fieldlength) = LFA4py.wlfalecr(q_unit, field, fieldlength)  
        else:
          print('Pb extraction for ' + field + '   fieldtype = ' + fieldtype)
          exit () 
      except RuntimeError:
   #     print('read_lfa_field : ' + field + ' not found in file') 
        raise NameError('FieldNotFound')
        return
     
      if extract: 
         if fieldlength == self.nlev*self.nbpoint:
            ideb = (self.point-1)*self.nlev
            ifin = self.point*self.nlev
         elif fieldlength == self.nlevp1*self.nbpoint:
            ideb = (self.point-1)*self.nlevp1
            ifin = self.point*self.nlevp1
         elif fieldlength == self.nbpoint:
            ideb = self.point-1
            ifin = self.point
         else:
            print('Pb extraction for ' + field + '   fieldlength = ',fieldlength)
            exit()           
         return data[ideb:ifin]
      else:
         return data   
 
#-----------------------------------------------------------------------------------#        
   def position(self,q_unit,ipoint):

      doc = 'DOCD' + '{:03d}'.format(ipoint)
      try:
        docd = self.read_lfa_field(q_unit,doc)
      except  NameError:
        if self.debug : print('Routine position ... ' + doc + 'not found in file')
        return

      lon_rad = docd[2]
      sin_lat = docd[3]
      lat_rad = np.arcsin(sin_lat)
      cos_lat = np.cos(lat_rad)
      lonr_rad = docd[6]
      sin_latr = docd[7]
      latr_rad = np.arcsin(sin_latr)
      cos_latr = np.cos(latr_rad)

      dist = self.RA*np.arccos(sin_lat*sin_latr + cos_lat*cos_latr*np.cos(lon_rad-lonr_rad))
      lon = lon_rad*180./np.pi 
      if lon>180. : lon = lon - 360.
      lat = lat_rad*180./np.pi
      lonr = lonr_rad*180./np.pi
      latr = latr_rad*180./np.pi
      return(lon,lat,lonr,latr,dist)
