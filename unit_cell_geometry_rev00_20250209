#
# TITLE: unit cell geometry
# AUTHOR: Hyunseung Yoo
# PURPOSE:
# REVISION:
# REFERENCE:
#

import numpy as np
import matplotlib.pyplot as plt


#
# CLASS: UNIT_CELL
#

class UNIT_CELL:

    def __init__(self, uc_name, ang_div, cut_cd, ch_recess):
        # user input 1
        self.uc_name = uc_name          # string

        # user input 2
        self.ang_div = ang_div          # number of division
        self.d_ang = 2.0*np.pi / self.ang_div
        self.gt_ang = np.arange(0.0, 2.0*np.pi, self.d_ang)

        # user input 3
        self.cut_cd_half = cut_cd / 2.0     # nm

        # user input 4
        self.ch_recess = ch_recess         # angstrom


    def ellipse_plug(self, gt_x_cd, gt_y_cd):
        # user input
        self.gt_x_cd = gt_x_cd      # nm
        self.gt_y_cd = gt_y_cd      # nm

        # calculations 1
        self.gt_r_x = (self.gt_x_cd / 2.0)     #nm
        self.gt_r_y = (self.gt_y_cd / 2.0)     # nm
        self.gt_x = self.gt_r_x * np.cos(self.gt_ang)      # nm
        self.gt_y = self.gt_r_y * np.sin(self.gt_ang)      # nm
        self.gt_r = np.sqrt(self.gt_x*self.gt_x + self.gt_y*self.gt_y)

        # calculations 2
        self.gt_peri = np.sum(self.gt_r*self.d_ang)         # nm


    def ponoa_alo(self, alo_thk):
        # user input
        self.alo_thk = alo_thk      # angstrom

        # calculations 1
        self.alo_r_x = (self.gt_r_x - self.alo_thk /10.0)       # nm
        self.alo_r_y = (self.gt_r_y - self.alo_thk /10.0)       # nm
        self.alo_x = self.alo_r_x * np.cos(self.gt_ang)         # nm
        self.alo_y = self.alo_r_y * np.sin(self.gt_ang)         # nm
        self.alo_r = np.sqrt(self.alo_x*self.alo_x + self.alo_y*self.alo_y)     # nm

        # calculations 2
        self.alo_peri = np.sum(self.alo_r * self.d_ang)         # nm


    def ponoa_box(self, box_thk):
        # user input
        self.box_thk = box_thk      # angstrom

        # calculations 1
        self.box_r_x = (self.alo_r_x - self.box_thk /10.0)      # nm
        self.box_r_y = (self.alo_r_y - self.box_thk /10.0)      # nm
        self.box_x = self.box_r_x * np.cos(self.gt_ang)         # nm
        self.box_y = self.box_r_y * np.sin(self.gt_ang)         # nm
        self.box_r = np.sqrt(self.box_x * self.box_x + self.box_y * self.box_y)  # nm

        # calculations 2
        self.box_peri = np.sum(self.box_r * self.d_ang)         # nm


    def ponoa_ctn(self, ctn_thk):
        # user input
        self.ctn_thk = ctn_thk      # angstrom

        # calculations 1
        self.ctn_r_x = (self.box_r_x - self.ctn_thk /10.0)      # nm
        self.ctn_r_y = (self.box_r_y - self.ctn_thk /10.0)      # nm
        self.ctn_x = self.ctn_r_x * np.cos(self.gt_ang)         # nm
        self.ctn_y = self.ctn_r_y * np.sin(self.gt_ang)         # nm
        self.ctn_r = np.sqrt(self.ctn_x * self.ctn_x + self.ctn_y * self.ctn_y)  # nm

        # calculations 2
        self.ctn_x_cut_mask_pos = self.ctn_x > +self.cut_cd_half        # mask, positive x axis
        self.ctn_x_cut_mask_neg = self.ctn_x < -self.cut_cd_half        # mask, negative x axis

        self.ctn_x_pos = self.ctn_x * self.ctn_x_cut_mask_pos           # positive x axis
        self.ctn_x_neg = self.ctn_x * self.ctn_x_cut_mask_neg           # negative x axis

        self.ctn_y_pos = self.ctn_y * self.ctn_x_cut_mask_pos           # positive x axis
        self.ctn_y_neg = self.ctn_y * self.ctn_x_cut_mask_neg           # negative x axis

        self.ctn_x_cut = self.ctn_x_pos + self.ctn_x_neg
        self.ctn_y_cut = self.ctn_y_pos + self.ctn_y_neg

        self.ctn_r_cut = np.sqrt(self.ctn_x_cut * self.ctn_x_cut + self.ctn_y_cut * self.ctn_y_cut)  # nm

        # calculations 3
        self.ctn_peri = np.sum(self.ctn_r * self.d_ang)                 # nm
        self.ctn_peri_cut = np.sum(self.ctn_r_cut * self.d_ang)         # nm


    def ponoa_tox(self, tox_thk):
        # user input
        self.tox_thk = tox_thk      # angstrom

        # calculations 1
        self.tox_r_x = (self.ctn_r_x - self.tox_thk /10.0)      # nm
        self.tox_r_y = (self.ctn_r_y - self.tox_thk /10.0)      # nm
        self.tox_x = self.tox_r_x * np.cos(self.gt_ang)         # nm
        self.tox_y = self.tox_r_y * np.sin(self.gt_ang)         # nm
        self.tox_r = np.sqrt(self.tox_x * self.tox_x + self.tox_y * self.tox_y)  # nm

        # calculations 2
        self.tox_x_cut_mask_pos = self.tox_x > +self.cut_cd_half        # mask, positive x axis
        self.tox_x_cut_mask_neg = self.tox_x < -self.cut_cd_half        # mask, negative x axis

        self.tox_x_pos = self.tox_x * self.tox_x_cut_mask_pos           # positive x axis
        self.tox_x_neg = self.tox_x * self.tox_x_cut_mask_neg           # negative x axis

        self.tox_y_pos = self.tox_y * self.tox_x_cut_mask_pos           # positive x axis
        self.tox_y_neg = self.tox_y * self.tox_x_cut_mask_neg           # negative x axis

        self.tox_x_cut = self.tox_x_pos + self.tox_x_neg
        self.tox_y_cut = self.tox_y_pos + self.tox_y_neg

        self.tox_r_cut = np.sqrt(self.tox_x_cut * self.tox_x_cut + self.tox_y_cut * self.tox_y_cut)  # nm

        # calculations 3
        self.tox_peri = np.sum(self.tox_r * self.d_ang)                 # nm
        self.tox_peri_cut = np.sum(self.tox_r_cut * self.d_ang)         # nm


    def ponoa_ch(self, ch_thk):
        # user input
        self.ch_thk = ch_thk      # angstrom

        # calculations 1
        self.ch_r_x = (self.tox_r_x - self.ch_thk /10.0)      # nm
        self.ch_r_y = (self.tox_r_y - self.ch_thk /10.0)      # nm
        self.ch_x = self.ch_r_x * np.cos(self.gt_ang)         # nm
        self.ch_y = self.ch_r_y * np.sin(self.gt_ang)         # nm
        self.ch_r = np.sqrt(self.ch_x * self.ch_x + self.ch_y * self.ch_y)  # nm

        # calculations 2
        self.ch_x_cut_mask_pos = self.ch_x > +self.cut_cd_half + self.ch_recess /10.0   # mask, positive x axis
        self.ch_x_cut_mask_neg = self.ch_x < -self.cut_cd_half - self.ch_recess /10.0   # mask, negative x axis

        self.ch_x_pos = self.ch_x * self.ch_x_cut_mask_pos           # positive x axis
        self.ch_x_neg = self.ch_x * self.ch_x_cut_mask_neg           # negative x axis

        self.ch_y_pos = self.ch_y * self.ch_x_cut_mask_pos           # positive x axis
        self.ch_y_neg = self.ch_y * self.ch_x_cut_mask_neg           # negative x axis

        self.ch_x_cut = self.ch_x_pos + self.ch_x_neg
        self.ch_y_cut = self.ch_y_pos + self.ch_y_neg

        self.ch_r_cut = np.sqrt(self.ch_x_cut * self.ch_x_cut + self.ch_y_cut * self.ch_y_cut)  # nm

        # calculations 3
        self.ch_peri = np.sum(self.ch_r * self.d_ang)                 # nm
        self.ch_peri_cut = np.sum(self.ch_r_cut * self.d_ang)         # nm


    def calculation_results(self):
        print('Unit Cell Name = %s' % self.uc_name)
        print('PLUG CD (major / minor) = %.1f / %.1f nm' % (self.gt_x_cd, self.gt_y_cd))
        print('CUT CD = %.1f nm' % (self.cut_cd_half * 2.0))
        print('CH recess = %.1f Angstrom' % (self.ch_recess))
        print('PONOA thickness = %.1f/%.1f/%.1f/%.1f/%.1f Angstrom' % \
              (self.ch_thk, self.tox_thk, self.ctn_thk, self.box_thk, self.alo_thk))
        print('')
        print('gate-alo interface perimeter = %.1f nm' % self.gt_peri)
        print('alo-box  interface perimeter = %.1f nm' % self.alo_peri)
        print('box-ctn  interface perimeter = %.1f nm' % self.box_peri)
        print('ctn-tox  interface perimeter = %.1f nm (all)' % self.ctn_peri)
        print('ctn-tox  interface perimeter = %.1f nm (cut /side -> %.2f percent)' % \
              (self.ctn_peri_cut/2, self.ctn_peri_cut/2/self.ctn_peri*100.0))
        print('tox-ch   interface perimeter = %.1f nm (all)' % self.tox_peri)
        print('tox-ch   interface perimeter = %.1f nm (cut /side -> %.2f percent)' % \
              (self.tox_peri_cut/2, self.tox_peri_cut/2/self.tox_peri*100.0))
        print('ch-fill  interface perimeter = %.1f nm (all)' % self.ch_peri)
        print('ch-fill  interface perimeter = %.1f nm (cut /side -> %.2f percent)' % \
              (self.ch_peri_cut/2, self.ch_peri_cut/2/self.ch_peri*100.0))



#
# main
#

uc_pdb = UNIT_CELL(uc_name='PDB', ang_div=360.0, cut_cd=60.0, ch_recess=0.0)       # string, ea, nm, angstrom
uc_pdb.ellipse_plug(gt_x_cd=180.0, gt_y_cd=110.0)                   # nm
uc_pdb.ponoa_alo(alo_thk=23.0)                                      # angstrom
uc_pdb.ponoa_box(box_thk=70.0)                                      # angstrom
uc_pdb.ponoa_ctn(ctn_thk=54.0)                                      # angstrom
uc_pdb.ponoa_tox(tox_thk=48.0)                                      # angstrom
uc_pdb.ponoa_ch(ch_thk=70.0)                                      # angstrom
uc_pdb.calculation_results()


# visualization
fig, ax = plt.subplots(1, 1)

x = uc_pdb.gt_x
y = uc_pdb.gt_y
ax.plot(x, y, linewidth=3.0)

x = uc_pdb.alo_x
y = uc_pdb.alo_y
ax.plot(x, y, linewidth=2.0)

x = uc_pdb.box_x
y = uc_pdb.box_y
ax.plot(x, y, linewidth=2.0)

x = uc_pdb.ctn_x
y = uc_pdb.ctn_y
ax.plot(x, y, ':',  linewidth=1.0)

x = uc_pdb.ctn_x_cut
y = uc_pdb.ctn_y_cut
ax.plot(x, y, linewidth=2.0)

x = uc_pdb.tox_x
y = uc_pdb.tox_y
ax.plot(x, y, ':',  linewidth=1.0)

x = uc_pdb.tox_x_cut
y = uc_pdb.tox_y_cut
ax.plot(x, y, linewidth=2.0)

x = uc_pdb.ch_x
y = uc_pdb.ch_y
ax.plot(x, y, ':', linewidth=1.0)

x = uc_pdb.ch_x_cut
y = uc_pdb.ch_y_cut
ax.plot(x, y, linewidth=2.0)


ax.axis('equal')
ax.set_xlim(-150,150)
ax.set_ylim(-100,100)
ax.grid(ls=':')

plt.show()
