""" 
    TLA region for a single sample:
 
"""

# %% Imports
import os
import sys
import psutil
import gc
import math
import time
import torch

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from PIL import Image
from argparse import ArgumentParser

from tla_functions import mkdirs, filexists

if torch.cuda.is_available():
    ISCUDA = True
else:
    ISCUDA = False

Image.MAX_IMAGE_PIXELS = 600000000

MINAREA = 100  # pixel^2

#IMGEXT = '.png'
IMGEXT = '.pdf'

__version__ = "2.0.2"


# %% Private classes

class Progress:
    """Progress class"""
    
    def __init__(self, name, maxsteps=10, sub=False):
        
        self.label = name
        self.start = time.time()
        self.mbuse = [np.nan] * maxsteps
        self.step = 0
        self.mbuse[0] = memuse()
        
        if sub:
            self.arrow = "====> " 
        else:
            self.arrow = "==> " 
            
            
    def dostep(self, display, lab=''):
        
        self.runtime = time.time() - self.start
        self.step = self.step + 1
        self.mbuse[self.step] = memuse()

        if display:
            print(self.arrow + self.label + \
                  "; STEP: " + str(self.step) + "; " + lab +'\n'\
                  " - Run Time: " + ftime(self.runtime) + "[HH:MM:SS]" +\
                  " - Mem used: " + str(self.mbuse[self.step] ) + "[MB]")
            

class Study:
    """Study class"""

    def __init__(self, study, pth):

        from tla_functions import tofloat, toint

        # loads arguments for this study
        self.name = study['name']
        self.raw_path = os.path.join(pth, study['raw_path'])

        # loads samples table for this study
        f = filexists(os.path.join(self.raw_path, study['raw_samples_table']))
        
        self.samples = pd.read_csv(f)

        # sets path for processed data
        self.dat_path = mkdirs(os.path.join(pth, study['data_path']))

        # list of processed samples
        self.done_list = os.path.join(self.dat_path, 'setup_done_samples.csv')

        # scale parameters
        self.factor = 1.0
        if 'factor' in study:
            self.factor = study['factor']
        self.scale = tofloat(study['scale']/self.factor)
        self.units = study['units']
        self.reg_factor = 4.0
        if 'reg_factor' in study:
            self.reg_factor = study['reg_factor']
            
        if 'image_file' in self.samples:
            aux = self.samples.image_file.astype(str)
            self.samples.image_file = aux.str.replace('nan', '', regex=False)
        else:
            self.samples.image_file = ''
            
        if 'mask_file' in study:
            aux = self.samples.mask_file.astype(str)
            self.samples.image_file = aux.str.replace('nan', '', regex=False)
        else:
            self.samples.mask_file = ''

        # the size of quadrats (to pixels as multiple of 10) (long scale)
        aux = np.rint(10*np.ceil((study['binsiz']/self.scale)/10))
        self.binsiz = toint(aux)
        
        # bandwidth size for convolutions (short scale)
        self.bw = toint(np.rint(self.binsiz/20))

        # loads classes info
        f = filexists(os.path.join(self.raw_path, 
                                   study['raw_reg_classes_table']))
        classes = pd.read_csv(f)
        classes['class'] = classes['class'].astype(str)
        self.classes = classes

        # creates tables for output
        self.allstats_out = pd.DataFrame()
        self.allpops_out = pd.DataFrame()


class Slide:
    """Slide class."""

    def __init__(self, i, study):

        # creates sample object
        self.tbl = study.samples.iloc[i].copy()  # table of parameters
        self.sid = self.tbl.sample_ID            # sample ID
        self.classes = study.classes             # region classes
        self.msg = "==> Slide [" + str(i + 1) + \
              "/" + str(len(study.samples.index)) + \
              "] :: Slide ID <- " + self.sid + " >>> pre-processing..."
        
        # paths
        self.raw_path = study.raw_path
        self.dat_path = study.dat_path

        # raw data file
        f = os.path.join(self.raw_path, self.tbl.coord_file)
        self.raw_cell_data_file = filexists(f)
        
        # roi mask image
        self.roi_file = os.path.join(self.raw_path, self.tbl.roi_file)
        # makes sure the path for roi file exists
        _ = mkdirs(os.path.dirname(self.roi_file))
        
        # reg mask image
        self.reg_file = filexists(os.path.join(self.raw_path, 
                                               self.tbl.region_file))
        
        # the size of quadrats and subquadrats
        self.binsiz = study.binsiz
        
        # bandwidth size for convolutions
        self.bw = study.bw

        # scale parameters
        self.factor = study.factor
        self.scale = study.scale
        self.units = study.units
        self.reg_factor = study.reg_factor

        # other sample attributes (to be filled later)
        self.imshape = []                # shape of landscape
        
        # if ROI should be split 
        self.split_rois = self.tbl.split_rois
        
        # sample arrays
        self.img = []                    # image array (slide image)
        self.msk = []                    # mask array  (blobs mask)
        self.roiarr = []                 # ROI array   (roi mask array)
        self.reg = []                    # region labels array
        
        # slide image
        if (self.tbl.image_file == ''):
            self.raw_imfile = ''
            self.isimg = False
        else:
            aux = os.path.join(study.raw_path, self.tbl.image_file)
            self.raw_imfile = aux
            self.isimg = True

        # blob mask image
        if (self.tbl.mask_file == ''):
            self.raw_mkfile = ''
            self.ismk = False
        else:
            aux = os.path.join(study.raw_path, self.tbl.mask_file)
            self.raw_mkfile = aux
            self.ismk = True
        
        # list of section labels in ROI array
        self.roi_sections = []


    def setup_data(self, edge=0):
        """
        Load region data and convert array to same dimensions as coordinate 
        data, also splits into sub-samples according to ROIs

        Args:
            - edge (int, optional):size [pix] of edge around data extremes 
             in rasters. Defaults to 0.
             
        Returns:
            None.

        """
        
        from skimage.transform import rescale
        from tla_functions import toint, tobit
        
        # load region array and rescale it
        aux = np.load(self.reg_file).T
        auy = rescale(aux, self.reg_factor*self.factor, anti_aliasing=True)
        # set pixel values to same range as original array
        regin = toint(np.round(auy*np.nanmax(aux)/np.nanmax(auy)))

        # read coordinate data
        cxy = pd.read_csv(self.raw_cell_data_file)[['class', 'x', 'y']]
        cxy['class'] = cxy['class'].astype(str)
        
        # updates coordinae values by conversion factor (from pix to piy)
        cxy.x = toint(cxy.x*self.factor)
        cxy.y = toint(cxy.y*self.factor)
        
        # gets extreme pixel values
        xmin, xmax = toint(np.min(cxy.x)), toint(np.max(cxy.x))
        ymin, ymax = toint(np.min(cxy.y)), toint(np.max(cxy.y))
        
        # limits for image cropping
        rmin = toint(np.nanmax([0, ymin - edge]))
        cmin = toint(np.nanmax([0, xmin - edge]))
        rmax = ymax + edge
        cmax = xmax + edge

        dr = rmax - rmin
        dc = cmax - cmin
        if (np.isnan(dr) or np.isnan(dc) or (dr <= 0) or (dc <= 0)):
            print("ERROR: data file " + self.raw_cell_data_file +
                  " is an empty or invalid landscape!")
            sys.exit()
        imshape = [toint(dr + 1), toint(dc + 1)]

        # shifts coordinates
        cell_data = xyShift(cxy, imshape, [rmin, cmin], self.scale)

        # crop region raster
        rmax_ = np.min([rmax, regin.shape[0]])
        cmax_ = np.min([cmax, regin.shape[1]])
        reg = np.zeros((imshape[0], imshape[1]), dtype='uint8')
        reg[0:(rmax_ - rmin),
            0:(cmax_ - cmin)] = tobit(regin[rmin:rmax_, cmin:cmax_,])
        self.regarr = reg

        self.cell_data = cell_data.reset_index(drop=True)
        self.imshape = [imshape[0], imshape[1]]
        
        
    def roi_mask(self, redo):
        """Generate ROI."""
        from tla_functions import filterCells, toint

        fout = self.roi_file
        if not os.path.exists(fout):
            # NOTE: if roi mask file does not exist, assumes this is a single
            #       sample slide, and filters out small disconnected regions
            from tla_functions import kdeMask_rois
            # gets a mask for the region that has cells inside
            self.roiarr = kdeMask_rois(self.cell_data,
                                       self.imshape,
                                       self.bw,
                                       100000,
                                       split=False,
                                       cuda=ISCUDA)
            np.savez_compressed(fout, roi=self.roiarr)
        else:
            aux = np.load(fout)
            self.roiarr = toint(aux['roi'])
            
        # filter out cells outside of ROI
        # adds a 'mask' tag for cells in different regions
        self.cell_data = filterCells(self.cell_data, self.roiarr)

        # list of roi sections
        self.roi_sections = np.unique(self.roiarr[self.roiarr > 0]).tolist()


class Sample:
    """Sample class."""

    def __init__(self, i, slide):

        from tla_functions import mkdirs, toint
        
        roii = slide.roi_sections[i]
        
        # get limits for this roi and crop reg array
        msk = (slide.roiarr == roii)

        roi = np.zeros_like(slide.roiarr, dtype='bool')
        roi[msk] = True
        
        reg = np.zeros_like(slide.regarr, dtype='int')
        reg[msk] = slide.regarr[msk].copy()
        
        # cell classes
        classes = slide.classes.copy()
        
        # drop regions that are excluded from study
        for k in classes.loc[classes['drop'], 'class'].tolist(): 
            reg[reg==toint(k)] = 0
        self.regarr = reg
        
        self.roiarr = roi # ROI array

        # creates sample object
        self.sid = str(slide.sid) + "_roi-" + str(roii) # sample ID
        self.tbl = slide.tbl.copy()                     # table of parameters
        self.tbl.sample_ID = self.sid
        self.msg = "====> Sample [" + str(i + 1) + \
              "/" + str(len(slide.roi_sections)) + \
              "] :: SID <- " + self.sid + " >>> pre-processing..."
        
        # drops classes that are excluded from study
        classes.drop(classes.loc[classes['drop']].index, inplace=True)
        classes.drop(columns=['drop'], inplace=True)
        classes.reset_index(inplace=True, drop=True)
        self.classes = classes

        # the size of quadrats and subquadrats
        self.binsiz = slide.binsiz

        # bandwidth size for convolutions
        self.bw = slide.bw

        # scale parameters
        self.factor = slide.factor
        self.scale = slide.scale
        self.units = slide.units
        self.reg_factor = slide.reg_factor

        # other sample attributes
        aux = slide.cell_data.copy()
        df = aux.loc[aux['mask'] == roii]
        self.cell_data = df                   # cell coordinates
        self.imshape = slide.imshape.copy()   # shape of accepted image
        
        # ploting limits (physical)
        self.ext = []

        # creates results folder and add path to sample tbl
        f = mkdirs(os.path.join(slide.dat_path,
                                'results', 'samples', self.sid))
        self.tbl['results_dir'] = 'results/samples/' + self.sid + '/'
        self.res_pth = f
        
        # sample arrays
        self.img = slide.img.copy()            # image array
        self.msk = slide.msk.copy()            # mask array (blobs)
        
        # creates images folder and add path to sample tbl
        if (slide.isimg):
            f = mkdirs(os.path.join(slide.dat_path, 'images'))
            self.tbl['image_file'] = 'images/' + self.sid + '_img.jpg'
            self.imfile = os.path.join(f, self.sid + '_img.jpg')
            self.isimg = True
        else:
            self.imfile = ''
            self.isimg = False
            
        # creates raster folder and add path to df
        f = mkdirs(os.path.join(slide.dat_path, 'rasters', self.sid))
        if (slide.ismk):
            self.tbl['mask_file'] = 'rasters/' + self.sid + '/' + \
                self.sid + '_mask.npz'
            self.mask_file = os.path.join(f, self.sid + '_mask.npz')
            self.ismk = True
        else:
            self.mask_file = ''
            self.ismk = False
            
        # raster file names
        self.raster_folder = f
        self.roi_file = os.path.join(f, self.sid + '_roi.npz')
        self.reg_file = os.path.join(f, self.sid + '_reg.npz')
        
        # classes info file
        f = os.path.join(self.res_pth, self.sid + '_reg_classes.csv')
        classes.to_csv(f, index=False)
        
        # output table
        self.out_file = os.path.join(self.res_pth, 
                                     self.sid + '_samples_reg.csv')


    def setup_data(self, edge):
        """
        Args:
            edge (int): size [pix] of edge around data extremes in rasters

        Returns:
            None.
        
        """
        
        from skimage import io
        from tla_functions import toint, tobit

        # resets coordinate values for cropped out sample
        cxy = self.cell_data[['class', 'col', 'row']]
        cxy.columns = ['class', 'x', 'y']
        
        # gets extreme pixel values
        xmin, xmax = toint(np.min(cxy.x)), toint(np.max(cxy.x))
        ymin, ymax = toint(np.min(cxy.y)), toint(np.max(cxy.y))

        # limits for image cropping
        rmin = toint(np.nanmax([0, ymin - edge]))
        cmin = toint(np.nanmax([0, xmin - edge]))
        rmax = ymax + edge
        cmax = xmax + edge

        imshape = [toint(rmax - rmin), toint(cmax - cmin)]

        # shifts coordinates
        cell_data = xyShift(cxy, imshape, [rmin, cmin], self.scale)

        # create croped versions of image and mask raster
        if self.isimg:
            ims = self.img.copy()
            rmax_ = np.min([rmax, ims.shape[0]])
            cmax_ = np.min([cmax, ims.shape[1]])
            img = np.zeros((imshape[0], imshape[1], 3), dtype='uint8')
            img[0:(rmax - rmin),
                0:(cmax - cmin), :] = tobit(ims[rmin:rmax_, cmin:cmax_, :])
            self.img = img
            io.imsave(self.imfile, self.img, check_contrast=False)
        else:
            self.img = []

        if self.ismk:
            msc = self.msk.copy()
            rmax_ = np.min([rmax, msc.shape[0]])
            cmax_ = np.min([cmax, msc.shape[1]])
            msk = np.zeros(imshape, dtype='uint8')
            msk[0:(rmax_ - rmin),
                0:(cmax_ - cmin)] = tobit(msc[rmin:rmax_, cmin:cmax_])
            [cell_data, msk_img] = getBlobs(cell_data, msk)
            np.savez_compressed(self.mask_file, roi=msk_img)
            self.msk = (msk_img > 0).astype('bool')
        else:
            self.msk = []

        aux = self.roiarr[rmin:rmax, cmin:cmax].copy()
        self.roiarr = aux
        fout = self.roi_file
        np.savez_compressed(fout, roi=self.roiarr)
        
        self.regarr = self.regarr[rmin:rmax, cmin:cmax]
        fout = self.reg_file
        np.savez_compressed(fout, reg=self.regarr)
        
        self.cell_data = cell_data.reset_index(drop=True)
        self.imshape = [toint(imshape[0]), toint(imshape[1])]
        
        # ploting limits (physical)
        self.ext = [0, np.round(self.imshape[1]*self.scale, 2),
                    0, np.round(self.imshape[0]*self.scale, 2)]


    def plotRegLandscape(self, out_pth):
        """
        Plot region landscape
  
        """
              
        from tla_functions import plotRGB, plotEdges, tofloat
        from matplotlib import colors
        import warnings
  
        nlevs = len(self.classes)
        raster = tofloat(self.regarr)
        raster[raster==0]= np.nan
        
        icticks = np.arange(nlevs) + 1
        cticks = self.classes['class_name'].tolist()
        ctitle = 'Region Categories'
        cmap = colors.ListedColormap(self.classes['class_color'].tolist())
        
        [ar, redges, cedges,  xedges, yedges] = plotEdges(self.imshape, 
                                                          4*self.binsiz, 
                                                          self.scale)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            
            # plots sample image
            fig, ax = plt.subplots(1, 1,
                                   figsize=(10*1, 0.5 + math.ceil(10*1/ar)),
                                   facecolor='w', edgecolor='k')
          
            im = plotRGB(ax, raster, self.units,
                         cedges, redges, xedges, yedges, fontsiz=18,
                         vmin=0.5, vmax=(nlevs + 0.5), cmap=cmap)
            #ax.grid(False)
            
            cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_ticks(icticks)
            cbar.set_ticklabels(cticks)
            cbar.ax.tick_params(labelsize=14)
            cbar.set_label(ctitle, size = 16, rotation=90, labelpad=1)
            ax.set_title('Tissue Regions', fontsize=24, y=1.02)
            #fig.subplots_adjust(hspace=0.4)
            fig.suptitle('Sample ID: ' + str(self.sid), fontsize=36, y=.95)
            #plt.tight_layout()
            fig.savefig(os.path.join(out_pth,
                                     self.sid +'_reg_landscape' + IMGEXT),
                        bbox_inches='tight', dpi=300)
            plt.close()
            
            # plots sample image detail
            # detail points
            lims2 = [2440, 3640, 2240, 3640]
            fig, ax = plt.subplots(1, 1,
                                   figsize=(10*1, 0.5 + math.ceil(10*1/ar)),
                                   facecolor='w', edgecolor='k')
          
            im = plotRGB(ax, raster, self.units,
                         cedges, redges, xedges, yedges, fontsiz=18,
                         vmin=0.5, vmax=(nlevs + 0.5), cmap=cmap)
            ax.set_xlim([lims2[0], lims2[1]])
            ax.set_ylim([lims2[2], lims2[3]])
            #ax.grid(False)
            cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_ticks(icticks)
            cbar.set_ticklabels(cticks)
            cbar.ax.tick_params(labelsize=14)
            cbar.set_label(ctitle, size = 16, rotation=90, labelpad=1)
            ax.set_title('Tissue Regions', fontsize=24, y=1.02)
            #fig.subplots_adjust(hspace=0.4)
            fig.suptitle('Sample ID: ' + str(self.sid), fontsize=36, y=.95)
            #plt.tight_layout()
            fig.savefig(os.path.join(out_pth,
                                     self.sid +'_reg_landscape_detail' + 
                                     IMGEXT),
                        bbox_inches='tight', dpi=300)
            plt.close()
            
            
            
  
            # the histogram of reg classes frequencies
            fig, ax = plt.subplots(1, 1, figsize=(5, 5),
                                   facecolor='w', edgecolor='k')
            ilabels, counts = np.unique(raster[~np.isnan(raster)],
                                        return_counts=True)
            ax.bar(ilabels, counts, align='center',
                   # alpha=0.5,
                   color=cmap(ilabels.astype(int)), edgecolor='k',)
            ax.set_title('Tissue Regions', fontsize=16, y=1.02)
            ax.set_xlabel(ctitle)
            ax.set_ylabel('Frequency')
            ax.set_xlim([0.5, (nlevs + 0.5)])
            ax.set_xticks(icticks)
            ax.set_xticklabels(cticks, rotation=90)
            ax.set_yscale('log')
            plt.tight_layout()
            fig.savefig(os.path.join(out_pth, 
                                     self.sid +'_reg_distribution' + IMGEXT),
                        bbox_inches='tight', dpi=300)
            plt.close()
            
    
    def landscapeAnalysis(self, minarea, redo, do_plots=False):
        
        samplcsv = os.path.join(self.res_pth, self.sid +'_reg_tbl.csv')
        
        if (redo or not os.path.exists(samplcsv)):

            import pylandstats as pls
            from skimage.morphology import convex_hull_image
            from skimage.measure import regionprops_table

            reg_pth = mkdirs(os.path.join(self.res_pth, 'regs'))

            self.regarr = cleanSmallRegions(self.regarr, minarea)
            
            # create pylandscape object and assing to landscape
            self.plsobj = pls.Landscape(self.regarr, res=(1, 1), 
                                        nodata=0, neighborhood_rule='8')
      
            # get adjacency matrix (as related to reg classes)
            adj_pairs = regAdjacency(self.plsobj)
            suff = '_reg_adjacency_odds.csv'
            adj_pairs.to_csv(os.path.join(reg_pth, self.sid + suff),
                             sep=',', index=False, header=True)

            # get all patch metrics 
            patch_metrics = getPatchMetrics(self.plsobj)
            suff = '_reg_patch_metrics.csv'
            patch_metrics.to_csv(os.path.join(reg_pth, self.sid + suff),
                                 sep=',', index=False, header=True)
            
            # get all class metrics 
            class_metrics = getClassMetrics(self.plsobj)
            
            class_metrics['axis_major_length'] = np.nan
            class_metrics['axis_minor_length'] = np.nan
            class_metrics['feret_diameter'] = np.nan
            class_metrics['convex_hull_area'] = np.nan
           
            prop_lst = ('axis_major_length', 
                        'axis_minor_length',
                        'area',
                        'feret_diameter_max')
            for sc in self.classes['class'].tolist():
                c = int(sc)
                chull = 1*convex_hull_image(1*(self.regarr == c))
                props = regionprops_table(chull,
                                          properties=prop_lst)
                
                val = props['axis_major_length'][0]
                class_metrics.loc[class_metrics['reg'] == c, 
                                  'axis_major_length'] = val
                val = props['axis_minor_length'][0]
                class_metrics.loc[class_metrics['reg'] == c, 
                                  'axis_minor_length'] = val
                val = props['area'][0]
                class_metrics.loc[class_metrics['reg'] == c, 
                                  'convex_hull_area'] = val
                val = props['feret_diameter_max'][0]
                class_metrics.loc[class_metrics['reg'] == c, 
                                  'feret_diameter'] = val
            
            suff = '_reg_class_metrics.csv'
            class_metrics.to_csv(os.path.join(reg_pth, self.sid + suff),
                                  sep=',', index=False, header=True)
            
            # get all landscape-level metrics 
            landscape_metrics = getLandscapeMetrics(self.plsobj)
            suff= '_reg_landscape_metrics.csv'
            landscape_metrics.to_csv(os.path.join(reg_pth, self.sid + suff),
                                     sep=',', index=True, 
                                     index_label="metric")
            
            if do_plots:
                # adjacency matrix heatmap
                adjacencyHeatmap(self.sid, adj_pairs, reg_pth)
                # plot patch metrics histograms
                _ = plotPatchHists(self.sid, patch_metrics, reg_pth)
                # plot class metrics histograms
                _ = plotClassHists(self.sid, class_metrics, reg_pth)
                
                # plot LME landscape
                self.plotRegLandscape(reg_pth)
       
            # get sample stats
            sample_out = getSampleStats(self.tbl, self.plsobj, adj_pairs,
                                        class_metrics, landscape_metrics)
            sample_out.to_csv(samplcsv, sep=',', index=False, header=True)
  

    def output(self, study, trun, memmax):
        """Record output."""
        samples_out = self.tbl.to_frame().T
        samples_out = samples_out.astype(study.samples.dtypes)
        samples_out['setup_runtime'] = trun
        samples_out['setup_maxMBs'] = memmax
        samples_out['index'] = np.nan
        samples_out.to_csv(self.out_file, index=False, header=True)

        # clean memory
        gc.collect()


# Private Functions

def memuse():
    """Memory use"""
    gc.collect()
    m = psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2
    return(round(m, 2))


def ftime(t):
    """Formated time string"""
    return(time.strftime('%H:%M:%S', time.gmtime(t)))
   

def xyShift(data, shape, ref, scale):
    """
    Shift coordinates and transforms into physical units.

    Args:
        - data (pandas): TLA dataframe of cell coordinates
        - shape (tuple): shape in pixels of TLA landscape
        - ref (tupe): reference location (upper-right corner)
        - scale (float): scale of physical units / pixel

    Returns:
        cell_data (pandas): data dable with shifted coordinates

    """
    from tla_functions import tofloat, toint

    cell_data = data.copy()

    # first drops duplicate entries
    # (shuffling first to keep a random copy of each duplicate)
    cell_data = cell_data.sample(frac=1.0).drop_duplicates(['x', 'y'])
    cell_data = cell_data.reset_index(drop=True)

    # generate a cell id
    # cell_data['cell_id'] = cell_data.index + 1

    # round pixel coordinates
    cell_data['col'] = toint(np.rint(cell_data['x']))
    cell_data['row'] = toint(np.rint(cell_data['y']))

    # shift coordinates to reference point
    cell_data['row'] = cell_data['row'] - ref[0]
    cell_data['col'] = cell_data['col'] - ref[1]

    # scale coordinates to physical units and transforms vertical axis
    cell_data['x'] = tofloat(cell_data['col']*scale)
    cell_data['y'] = tofloat((shape[0] - cell_data['row'])*scale)

    # drops data points outside of frames
    cell_data = cell_data.loc[(cell_data['row'] >= 0) &
                              (cell_data['row'] < shape[0]) &
                              (cell_data['col'] >= 0) &
                              (cell_data['col'] < shape[1])]

    return cell_data


def getBlobs(data, mask):
    """
    Get labels from blob regions mask, and assing them to the cell data.

    Args:
        data (pandas): TLA dataframe of cell coordinates
        mask (numpy): binary mask defining blob regions

    Returns:
        aux (pandas): dataframe with extra column for `blob`
        blobs_labels (numpy): array with blob labels

    """
    from tla_functions import toint
    from skimage import measure

    # get a binary image from the (thresholded) mask
    msk_img = np.zeros(mask.shape, dtype='uint8')
    msk_img[mask > 127] = 1

    # label blobs in mask image
    blobs_labels = measure.label(msk_img, background=0, connectivity=2)

    # get data coordinates and labels for mask
    rows, cols = np.where(msk_img > 0)
    msk_data = pd.DataFrame({'blob': blobs_labels[rows, cols],
                             'row': rows,
                             'col': cols})

    aux = pd.merge(data, msk_data, how="left",
                   on=["row", "col"]).fillna(0)
    aux['row'] = toint(aux['row'])
    aux['col'] = toint(aux['col'])
    aux['blob'] = toint(aux['blob'])

    blobs_labels = toint(blobs_labels)

    return (aux, blobs_labels)


# Space Statistics functions   

def regAdjacency(ls):
    """
    Gets reg points adjacency matrix

    Args:
        - ls: (pylandstats) landscape object
        - classes_file: (str) name of file with cell classes definitions
    """

    # get adjacency matrix (as related to reg classes)
    adj_pairs = ls.compute_total_adjacency_df()

    # number of (valid) adjacent pairs
    adj_pairs['total'] = adj_pairs.sum(axis=1)
    adj_pairs['reg'] = adj_pairs.index
    
    cols = list(adj_pairs)
    cols.insert(0, cols.pop(cols.index('reg')))
    adj_pairs = adj_pairs.loc[:, cols]
    adj_pairs.index.name = None
    
    adj_pairs['areax4'] = 0
    arr = ls.landscape_arr
    for i, c in enumerate(ls.classes):
        # odds per class
        adj_pairs.loc[adj_pairs['reg'] == c,
                      'areax4'] = int(4*np.sum(arr == int(c)))
        adj_pairs[c] = adj_pairs[c]/adj_pairs['total']

    return(adj_pairs)
 
    
def adjacencyHeatmap(sid, adj_pairs, res_pth ):
    
    import seaborn as sns
    
    mat = adj_pairs.copy()
    mat = mat.set_index('reg') 
    mat = mat[adj_pairs.reg]
    
    mat[mat==0] = 10**(np.floor(np.log10(np.nanmin(mat[mat>0]))))
    mat = np.log10(mat)
    sns.set(rc = {'figure.figsize':(10,10)})
    f = sns.heatmap(mat, linewidth=0.5, cmap="jet", vmax = 0,
                    cbar_kws={'label': r'$\log_{10}(F)$',
                              'shrink': 0.87,
                              'pad': 0.01},
                    xticklabels = 1, yticklabels = 1)
    f.set(title='Region adjacenty odds',
          xlabel=None, ylabel=None);
    f.invert_yaxis()
    f.set_aspect('equal')
    plt.yticks(rotation=0)
    f.figure.savefig(os.path.join(res_pth, 
                                  sid +'_reg_adjacency_odds' + IMGEXT),
                     bbox_inches='tight', dpi=300)
    plt.tight_layout()
    plt.close()

    return(0)


def cleanSmallRegions(arrin, minarea):
    
    from skimage.measure import label, regionprops_table
    
    arr = arrin.copy().astype(int)
    
    # label individual regions
    labls = label(arr, background=0, connectivity=2)
    
    # # This part uses a gaussian convolution to estimate density of regions
    # # Not a great solution but could allow us to merge many small and close
    # # together regions into larger ones that survive the size filtering
    # from tla_functions import gkern, fftconv2d
    # bw = np.sqrt(minarea/np.pi)
    # kern = gkern(bw)
    # # Compute the kernel density estimate
    # z1 = fftconv2d(arr==1, kern, cuda=False)
    # z2 = fftconv2d(arr==2, kern, cuda=False)
    # thr = 0.5
    # a1 = z1 > thr
    # l1 = label(a1, background=0, connectivity=2)
    # p1 = pd.DataFrame(regionprops_table(l1, properties=('label', 'area')))
    # m1 = np.isin(l1, p1['label'] )
    # a1[~m1] = 0
    # a2 = z2 > thr
    # l2 = label(a2, background=0, connectivity=2)
    # p2 = pd.DataFrame(regionprops_table(l2, properties=('label', 'area')))
    # m2 = np.isin(l2, p2['label'] )
    # a2[~m2] = 0
    # arr = 1*a1 + 2*a2
    # # make sure overlaps are taken as class=1
    # arr[a1>0] = 1
    
    # # This part uses a filter based on adjacencies of regions to delete
    # # small region by replacing them with the context around them
    # # This is VERY SLOW for samples that are very fragmented
    # from tla_functions import toint
    # import pylandstats as pls
    #
    # # label individual regions on filtered array
    # labls = label(arr, background=0, connectivity=2)
    # # gets pls object
    # ls = pls.Landscape(labls, res=(1, 1), nodata=0, neighborhood_rule='8')
    # # calculate adjacency mattrix for labeled regions
    # aux = ls.compute_total_adjacency_df()
    # aux = np.array(aux)
    # G = toint(aux > 0)
    # del aux
    # # filter out small regions: assigns a class value corresponding to a
    # # large neighboring region (picked at random)
    # p2 = props.loc[props['area'] < minarea]
    # for index, row in p2.iterrows():
    #     oldlab = int(row['label'])
    #     # gets all neighbors
    #     vals = (np.array(np.where(G[oldlab - 1, :] == 1)[0]) + 1).tolist()
    #     if len(vals) > 0:
    #         # finds area of all neighbors
    #         areas = props.loc[props['label'].isin(vals)]
    #         # pick a neighbor to replace the class value
    #         p3 = areas.loc[areas['area'] > minarea]
    #         if len(p3) > 0:
    #             newlab = np.random.choice(p3['label'].values)
    #             arr[labls==oldlab] = props.loc[props['label']==newlab, 'class']
    #         else:
    #             arr[labls==oldlab] = 0
    #     else:
    #         arr[labls==oldlab] = 0
    
    # calculate region properties
    props = regionprops_table(labls, properties=('label', 'area'))
    props = pd.DataFrame(props)
    psmall = props.loc[props['area'] < minarea]
    arr[np.isin(labls, psmall['label'])] = 0
    
    return(arr)

def getPatchMetrics(ls):

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')

        kws = {'area': {'hectares': False},
               'perimeter_area_ratio': {'hectares': False}}

        df = ls.compute_patch_metrics_df(metrics_kwargs=kws).reset_index()
        df['area_fraction'] = df['area']/ls.landscape_area
        df.rename(columns={'class_val': 'reg'}, inplace=True)
        
        # shift column 'LME' to first position
        df.insert(0, 'reg', df.pop('reg'))

    return(df)


def getClassMetrics(ls):
    ###########################################################################
    # Patch class-level metrics

    # Stadistics across patches for each class
    # `_mn` -> mean
    # `_am` -> area-weighted mean
    # `_md` -> median
    # `_ra` -> range
    # `_sd` -> standard deviation
    # `_cv` -> coefficient of variation

    # class_metrics.columns
    # ['total_area', 'proportion_of_landscape', 'number_of_patches',
    # 'patch_density', 'largest_patch_index',
    # 'total_edge', 'edge_density',
    # 'landscape_shape_index', 'effective_mesh_size',
    # 'area_mn', 'area_am', 'area_md', 'area_ra', 'area_sd', 'area_cv',
    # 'perimeter_mn', 'perimeter_am', 'perimeter_md', 'perimeter_ra',
    # 'perimeter_sd', 'perimeter_cv', 'perimeter_area_ratio_mn',
    # 'perimeter_area_ratio_am', 'perimeter_area_ratio_md',
    # 'perimeter_area_ratio_ra', 'perimeter_area_ratio_sd',
    # 'perimeter_area_ratio_cv', 'shape_index_mn', 'shape_index_am',
    # 'shape_index_md', 'shape_index_ra', 'shape_index_sd', 'shape_index_cv',
    # 'fractal_dimension_mn', 'fractal_dimension_am', 'fractal_dimension_md',
    # 'fractal_dimension_ra', 'fractal_dimension_sd','fractal_dimension_cv',
    # 'euclidean_nearest_neighbor_mn', 'euclidean_nearest_neighbor_am',
    # 'euclidean_nearest_neighbor_md', 'euclidean_nearest_neighbor_ra',
    # 'euclidean_nearest_neighbor_sd', 'euclidean_nearest_neighbor_cv']

    import warnings
    
    kws = {'total_area': {'hectares': False},
           'perimeter_area_ratio': {'hectares': False},
           'patch_density': {'hectares': False},
           'edge_density': {'hectares': False},
           'effective_mesh_size': {'hectares': False}}
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        df = ls.compute_class_metrics_df(metrics_kwargs=kws).reset_index()
        
        df['reg'] = df['class_val']
        
    # shift column 'LME' to first position
    df.insert(0, 'reg', df.pop('reg'))
        
    return(df)


def getLandscapeMetrics(ls):
    ###########################################################################
    # Landscape-level metrics: aggregating over all patches of the landscape
 
    # Stadistics across patches 
    # `_mn` -> mean
    # `_am` -> area-weighted mean
    # `_md` -> median
    # `_ra` -> range
    # `_sd` -> standard deviation
    # `_cv` -> coefficient of variation (=sigma/mu)
 
    # class_metrics.columns 
    #['total_area', 'number_of_patches', 'patch_density', 
    # 'largest_patch_index', 'total_edge', 'edge_density', 
    # 'landscape_shape_index', 'effective_mesh_size', 
    # 'area_mn', 'area_am', 'area_md', 'area_ra', 'area_sd', 'area_cv', 
    # 'perimeter_mn', 'perimeter_am', 'perimeter_md', 'perimeter_ra', 
    # 'perimeter_sd', 'perimeter_cv', 'perimeter_area_ratio_mn', 
    # 'perimeter_area_ratio_am', 'perimeter_area_ratio_md', 
    # 'perimeter_area_ratio_ra', 'perimeter_area_ratio_sd', 
    # 'perimeter_area_ratio_cv', 
    # 'shape_index_mn', 'shape_index_am', 'shape_index_md', 'shape_index_ra', 
    # 'shape_index_sd', 'shape_index_cv',
    # 'fractal_dimension_mn', 'fractal_dimension_am', 'fractal_dimension_md', 
    # 'fractal_dimension_ra', 'fractal_dimension_sd','fractal_dimension_cv', 
    # 'euclidean_nearest_neighbor_mn', 'euclidean_nearest_neighbor_am', 
    # 'euclidean_nearest_neighbor_md',
    # 'euclidean_nearest_neighbor_ra', 'euclidean_nearest_neighbor_sd', 
    # 'euclidean_nearest_neighbor_cv',
    # 'contagion', 'shannon_diversity_index']
    
    import warnings
    
    kws = {'total_area': {'hectares': False},
           'edge_density': {'hectares': False},
           'perimeter_area_ratio': {'hectares': False},
           'patch_density': {'hectares': False}}
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        df = ls.compute_landscape_metrics_df(metrics_kwargs=kws).T.squeeze()

    return(df.rename('value'))


def getSampleStats(sample, ls, adj_pairs, class_metrics, landscape_metrics):
    
    # gather sample statistics
    sample_out = sample.copy()
    
    npairs = np.sum(adj_pairs['total'].values)
    npairs_eq = 0
    for i, c in enumerate(adj_pairs['reg']):
        npairs_eq += (adj_pairs['total']*adj_pairs[c]).tolist()[i]
     
    idx = np.logical_or(sample_out.index.str.endswith('_file'),
                         sample_out.index.str.endswith('_dir'))
    sample_out = sample_out.drop(sample.index[idx].tolist())
    sample_out['num_reg_classes'] = len(ls.classes)
    sample_out['landscape_area'] = ls.landscape_area
    sample_out['adjacency_index'] = npairs_eq/npairs
    
    # Interspersion (or Contagion): measure of aggregation that 
    # represents the probability that two random adjacent cells belong 
    # to the same class.
    # => Approaches '0' when the classes are maximally disaggregated 
    #    (i.e., every cell is a patch of a different class) and interspersed 
    #    (i.e., equal proportions of all pairwise adjacencies)
    # => Approaches '100' when the whole landscape consists of a single patch.
    sample_out['contagion'] = landscape_metrics['contagion']

    # Shannon Diversity Index (Entropy): measure of diversity that reflects 
    # the number of classes present in the landscape as well as the relative 
    # abundance of each class.
    # => Approaches 0 when the entire landscape consists of a single patch.
    # => Increases as the number of classes increases and/or the proportional
    #    distribution of area among classes becomes more equitable.
    sample_out['Shannon'] = landscape_metrics['shannon_diversity_index']

    # Landscape Shape Index: measure of class aggregation that provides a 
    # standardized measure of edginess that adjusts for the size of the 
    # landscape.
    # => Equals 1 when the entire landscape consists of a single patch of the 
    #    corresponding class,
    # => Increases without limit as the patches become more disaggregated and 
    #    uniformly mixed.
    sample_out['shape_index'] = landscape_metrics['landscape_shape_index']
    
    # Simpson Diversity Index (1 - Dominance): measure of diversity that 
    # reflects the dominance of species, taking into account the relative 
    # abundance of each class.
    # => Equals 1 large diversity (more even distribution with no species 
    #    dominating the spectrum) 
    # => Equals 0 small diversity (the entire landscape consists of a single 
    #    patch of the corresponding class)
    ni = np.array(class_metrics.total_area)
    N  = np.sum(ni)
    sample_out['Simpson'] = np.sum(ni*(ni-1))/(N*(N-1))

    # Other metrics
    sample_out['num_patches'] = landscape_metrics['number_of_patches']
    sample_out['patch_density'] = landscape_metrics['patch_density']
    sample_out['total_edge'] = landscape_metrics['total_edge']
    sample_out['edge_density'] = landscape_metrics['edge_density']
    sample_out['largest_patch_ind'] = landscape_metrics['largest_patch_index']

    return(sample_out.to_frame().T)
     

def plotPatchHists(sid, df, res_pth):
    
    import matplotlib as mpl

    lmes = pd.unique(df['reg']).tolist()
    lab = ['Region class: ' + str(i) for i in lmes]

    # plot some patch metrics distributions
    fig, axs = plt.subplots(2, 2, figsize=(10, 10),
                            facecolor='w', edgecolor='k')

    col = [mpl.cm.jet(x) for x in np.linspace(0, 1, len(lmes))]

    for i, lme in enumerate(lmes):

        mif = np.log10(np.min(df.area_fraction))
        maf = np.log10(np.max(df.area_fraction))
        dat = df.loc[df['reg'] == lme].area_fraction.apply(np.log10)
        hist, bins = np.histogram(dat,
                                  bins=np.arange(mif,
                                                 maf,
                                                 (maf-mif)/21),
                                  density=False)
        x = (bins[1:] + bins[:-1])/2
        axs[0, 0].plot(x, hist, label=lab[i], color=col[i])

        mir = np.log10(np.min(df.perimeter_area_ratio))
        mar = np.log10(np.max(df.perimeter_area_ratio))
        dat = df.loc[df['reg'] == lme].perimeter_area_ratio.apply(np.log10)
        hist, bins = np.histogram(dat,
                                  bins=np.arange(mir,
                                                 mar,
                                                 (mar-mir)/21),
                                  density=False)
        x = (bins[1:] + bins[:-1])/2
        axs[0, 1].plot(x, hist, label=lab[i], color=col[i])

        mis = np.log10(np.min(df.shape_index))
        mas = np.log10(np.max(df.shape_index))
        dat = df.loc[df['reg'] == lme].shape_index.apply(np.log10)
        hist, bins = np.histogram(dat,
                                  bins=np.arange(mis,
                                                 mas,
                                                 (mas-mis)/21),
                                  density=False)

        x = (bins[1:] + bins[:-1]) / 2
        axs[1, 0].plot(x, hist, label=lab[i], color=col[i])

        mie = np.log10(np.min(df.euclidean_nearest_neighbor))
        mae = np.log10(np.max(df.euclidean_nearest_neighbor))
        aux = df.loc[df['reg'] == lme]
        dat = aux.euclidean_nearest_neighbor.apply(np.log10)
        hist, bins = np.histogram(dat,
                                  bins=np.arange(mie,
                                                 mae,
                                                 (mae-mie)/21),
                                  density=False)
        x = (bins[1:] + bins[:-1])/2
        axs[1, 1].plot(x, hist, label=lab[i], color=col[i])

    #axs[0, 0].set_xticks(np.arange(-10, 10, 1))
    axs[0, 0].set_xlim(mif, maf)
    axs[0, 0].locator_params(tight=True, nbins=5)
    axs[0, 0].get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(
        lambda x, p: r'$10^{%.1f}$' % x))
    axs[0, 0].set_xlabel('Area fraction of patches (log)')
    axs[0, 0].set_ylabel('Frequency')

    #axs[0, 1].set_xticks(np.arange(-10, 10, 1))
    axs[0, 1].set_xlim(mir, mar)
    axs[0, 1].locator_params(tight=True, nbins=5)
    axs[0, 1].get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(
        lambda x, p: r'$10^{%.1f}$' % x))
    axs[0, 1].set_xlabel('Perimeter-Area ratio of patches (log)')
    axs[0, 1].set_ylabel('Frequency')

    axs[1, 0].set_xlim(mis, mas)
    axs[1, 0].locator_params(tight=True, nbins=5)
    axs[1, 0].get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(
        lambda x, p: r'$10^{%.1f}$' % x))
    axs[1, 0].set_xlabel('Shape index of patches (log)')
    axs[1, 0].set_ylabel('Frequency')

    axs[1, 1].set_xlim(mie, mae)
    axs[1, 1].locator_params(tight=True, nbins=5)
    axs[1, 1].get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(
        lambda x, p: r'$10^{%.1f}$' % x))
    axs[1, 1].set_xlabel('Euclidean Nearest Neighbor of patches (log)')
    axs[1, 1].set_ylabel('Frequency')


    fig.legend(lmes, 
               loc='center right', 
               bbox_to_anchor=(1.05,0.5), 
               ncol=1, 
               bbox_transform=fig.transFigure)
    fig.subplots_adjust(hspace=0.4)
    fig.suptitle('Sample ID: ' + str(sid), fontsize=24, y=.95)
    #plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(os.path.join(res_pth, sid +'_reg_patch_metrics_hists' + 
                             IMGEXT),
                bbox_inches='tight', dpi=300)
    plt.close()

    return(fig)


def plotClassHists(sid, df, res_pth):
    # plot some patch metrics distributions
    
    import matplotlib as mpl

    lmes = pd.unique(df['reg']).tolist()
    cols = [mpl.cm.jet(x) for x in np.linspace(0, 1, len(lmes))]

    fig, axs = plt.subplots(2, 2, figsize=(10, 10),
                            facecolor='w', edgecolor='k')

    axs[0, 0].bar(lmes, df.patch_density, color=cols, align='center')
    axs[0, 0].set_title('Patch Density')
    axs[0, 0].set_xticks(lmes)
    axs[0, 0].set_xticklabels(lmes, rotation=90, fontsize=8)
    axs[0, 1].bar(lmes, df.largest_patch_index, color=cols, align='center')
    axs[0, 1].set_title('Largest Patch Index')
    axs[0, 1].set_xticks(lmes)
    axs[0, 1].set_xticklabels(lmes, rotation=90, fontsize=8)
    axs[1, 0].bar(lmes, df.edge_density, color=cols, align='center')
    axs[1, 0].set_title('Edge Density')
    axs[1, 0].set_xticks(lmes)
    axs[1, 0].set_xticklabels(lmes, rotation=90, fontsize=8)
    axs[1, 1].bar(lmes, df.landscape_shape_index, color=cols, align='center')
    axs[1, 1].set_title('Landscape Shape Index')
    axs[1, 1].set_xticks(lmes)
    axs[1, 1].set_xticklabels(lmes, rotation=90, fontsize=8)

    fig.subplots_adjust(hspace=0.4)
    fig.suptitle('Sample ID: ' + str(sid), fontsize=24, y=.95)
    fig.savefig(os.path.join(res_pth, sid +'_reg_class_metrics_hists' + 
                             IMGEXT),
                bbox_inches='tight', dpi=300)
    #plt.tight_layout()
    plt.close()

    # get class coverage graphs
    fig2, axs = plt.subplots(1, 2, figsize=(10, 5),
                             facecolor='w', edgecolor='k')

    y = df.total_area
    f = 100.*y/y.sum()
    patches, texts = axs[0].pie(y, colors=cols, startangle=10)
    axs[0].axis('equal')
    axs[0].set_title('Proportion of landscape')
    labs = ['{0} - {1:1.2f} %'.format(i, j) for i, j in zip(lmes, f)]
    axs[0].legend(patches, labs, 
                  loc='center left',
                  bbox_to_anchor=(0.0, 0.5), 
                  fontsize=8,
                  prop={'size': 8},
                  bbox_transform=fig2.transFigure)

    y = df.number_of_patches
    f = 100.*y/y.sum()
    patches, texts = axs[1].pie(y, colors=cols, startangle=0)
    axs[1].axis('equal')
    axs[1].set_title('Number of Patches')
    labs = ['{0} - {1:1.2f} %'.format(i, j) for i, j in zip(lmes, f)]
    axs[1].legend(patches, labs, 
                  loc='center right',
                  bbox_to_anchor=(1.0, 0.5), 
                  fontsize=8,
                  prop={'size': 8},
                  bbox_transform=fig2.transFigure)

    fig2.subplots_adjust(hspace=0.5)
    fig2.suptitle('Sample ID: ' + str(sid), fontsize=24, y=.95)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig2.savefig(os.path.join(res_pth, sid +'_reg_class_coverage' + IMGEXT),
                 bbox_inches='tight', dpi=300)
    plt.close()

    return([fig, fig2])

# %% Main function

def main(args):
    """TLA Setup Main."""

    # %% STEP 0: start, checks how the program was launched
    debug = False
    try:
        args
    except NameError:
        #  if not running from the CLI, run in debug mode
        debug = True

    if debug:
        # running from the IDE
        import seaborn as sns
        # path of directory containing this script
        main_pth = os.path.dirname(os.getcwd())
        f = os.path.join(main_pth, 'test_set_reg.csv')
        REDO = True
        GRPH = True
        CASE = 0   # case number to process

    else:
        # running from the CLI, (eg. using the bash script)
        # path to working directory (above /scripts)
        main_pth = os.getcwd()
        f = os.path.join(main_pth, args.argsfile)
        REDO = args.redo
        GRPH = args.graph
        CASE = args.casenum

    argsfile = filexists(f)

    print("-------------------------------------------------")    
    print("==> The working directory is: " + main_pth)
    print("==> Is CUDA available: " + str(ISCUDA))

    # NOTE: only the FIRST line in the argument table will be used
    study = Study(pd.read_csv(argsfile).iloc[0], main_pth)

    # STEP 1: create slide object
    slide = Slide(CASE, study)
    print(slide.msg)

    # tracks time and memory usage
    progress = Progress(slide.sid)
    progress.dostep(debug, 'Slide object created')


    # %% STEP 2: loads and format coordinate data for slide
    slide.setup_data()

    if debug:
        plt.close('all')
        plt.figure()
        sns.scatterplot(data=slide.cell_data, s=1,
                        x='x', y='y', hue='class')
        plt.legend(bbox_to_anchor=(1.02, 1), 
                   loc='upper left',
                   borderaxespad=0,
                   markerscale=5)
        plt.title("Slide ID: " + slide.sid)
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        del ax

        from matplotlib import colors
        plt.figure()
        cmap = colors.ListedColormap(['#ffffff'] + 
                                     slide.classes.class_color.tolist())
        bounds = np.arange(cmap.N + 1)
        norm = colors.BoundaryNorm(bounds, cmap.N)
        img = plt.imshow(slide.regarr, 
                         interpolation='nearest', 
                         origin='lower',
                         cmap=cmap, norm=norm)
        cbar = plt.colorbar(img, cmap=cmap, norm=norm, boundaries=bounds)
        tick_locs = bounds[:-1] + 0.5
        cbar.set_ticks(tick_locs)
        cbar.set_ticklabels(bounds[:-1])
        plt.title("Slide ID: " + slide.sid)
        plt.ylim(reversed(plt.ylim()))
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        del ax


    # STEP 3: load or calculate ROI mask for regions with cells. 
    #            ROIs are large unconnected sections of tissue 
    slide.roi_mask(REDO)

    if debug:
        plt.close('all')
        plt.figure()
        ext = [0, np.round(slide.imshape[1]*slide.scale, 2),
               0, np.round(slide.imshape[0]*slide.scale, 2)]
        plt.imshow(slide.roiarr, extent=ext)
        plt.title("Slide ID: " + slide.sid)
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        del ax

    progress.dostep(debug, 'ROIs loaded')

    # runtime up to here
    thead = progress.runtime
    print("==> ROI sections detected: " + str(slide.roi_sections))


    # STEP 4: for each ROI in the slide, create a new sample

    for i, roi in enumerate(slide.roi_sections):
        
        print("-------------------------------------------------")
        sample = Sample(i, slide)
        print(sample.msg)

        # tracks time and memory usage for this sample
        subprogress = Progress(sample.sid, sub=True)
        subprogress.dostep(debug, 'Sample object created')

        # if pre-processed files do not exist
        
        if (REDO or
                (not os.path.exists(sample.reg_file)) or
                (not os.path.exists(sample.roi_file)) or
                (not os.path.exists(sample.out_file))):

            #  STEP 5: setup coordinate data (pick out sample data)
            sample.setup_data(0)
            
            if debug:
                plt.close('all')
                plt.figure()
                sns.scatterplot(data=sample.cell_data, s=1,
                                x='x', y='y', hue='class')
                plt.legend(bbox_to_anchor=(1.02, 1), 
                           loc='upper left',
                           borderaxespad=0,
                           markerscale=5)
                plt.title("Sample_ID: " + sample.sid)
                ax = plt.gca()
                ax.set_aspect('equal', adjustable='box')
                del ax
                
                plt.figure()
                plt.imshow(sample.roiarr, extent=sample.ext)
                plt.title("ROI for Sample_ID: " + sample.sid)
                ax = plt.gca()
                ax.set_aspect('equal', adjustable='box')
                del ax
                
                from matplotlib import colors
                plt.figure()
                cols = ['#ffffff'] +  sample.classes.class_color.tolist()
                cmap = colors.ListedColormap(cols)
                bounds = np.arange(cmap.N + 1)
                norm = colors.BoundaryNorm(bounds, cmap.N)
                img = plt.imshow(sample.regarr, 
                                 interpolation='none', 
                                 origin='lower',
                                 cmap=cmap, norm=norm)
                cbar = plt.colorbar(img, cmap=cmap, 
                                    norm=norm, boundaries=bounds)
                tick_locs = bounds[:-1] + 0.5
                cbar.set_ticks(tick_locs)
                cbar.set_ticklabels(bounds[:-1])
                
                plt.title("Slide ID: " + sample.sid)
                plt.ylim(reversed(plt.ylim()))
                ax = plt.gca()
                ax.set_aspect('equal', adjustable='box')
                del ax
                
            subprogress.dostep(debug, 'Cell data formated')


            # STEP 6: pylandstats analysis
            # Regular metrics can be computed at the patch, class and
            # landscape level. For list of all implemented metrics see: 
            # https://pylandstats.readthedocs.io/en/latest/landscape.html
            
            sample.landscapeAnalysis(MINAREA, REDO, do_plots=GRPH)
            progress.dostep(debug, 'Landscape Analysis calculated')
            
            MINAREAum = np.round(MINAREA*sample.scale*sample.scale, 2)
            
            if debug:
                
                from matplotlib import colors
                plt.figure()
                cols = ['#ffffff'] +  sample.classes.class_color.tolist()
                cmap = colors.ListedColormap(cols)
                bounds = np.arange(cmap.N + 1)
                norm = colors.BoundaryNorm(bounds, cmap.N)
                img = plt.imshow(sample.regarr, 
                                 interpolation='none', 
                                 origin='lower',
                                 cmap=cmap, norm=norm)
                cbar = plt.colorbar(img, cmap=cmap, 
                                    norm=norm, boundaries=bounds)
                tick_locs = bounds[:-1] + 0.5
                cbar.set_ticks(tick_locs)
                cbar.set_ticklabels(bounds[:-1])
                plt.title("Slide ID: " + sample.sid + 
                          "\n(filtered @ " + str(MINAREA) + "[pix]^2 :: " +
                          str(MINAREAum) + 
                          sample.units + "^2)")
                plt.ylim(reversed(plt.ylim()))
                ax = plt.gca()
                ax.set_aspect('equal', adjustable='box')
                del ax
                
                #plt.savefig(os.path.join(sample.res_pth, 
                #             sample.sid + 'regions_small' + IMGEXT),
                #bbox_inches='tight', dpi=900)
            

        # LAST step: saves study stats results for sample
        memmax = np.max((np.nanmax(subprogress.mbuse), 
                         np.nanmax(progress.mbuse)))
        trun = ftime(thead + subprogress.runtime)
        print('====> Sample: ' + sample.sid + ' finished. Time elapsed: ', 
              trun, '[HH:MM:SS]')
        print("====> Max memory used: " + str(memmax) + "[MB]")
        
        sample.output(study, trun, memmax)
        
    # finish slide       
    with open(study.done_list, 'a') as f:
        f.write(slide.sid + '\n')
    
    #  no ROI        
    if (np.sum(slide.roiarr)==0):
        print("WARNING: Slide <" + slide.sid + "> doesn't have any valid ROIs")
        sys.exit()
        
    progress.dostep(False)
    trun = ftime(progress.runtime)
    print("-------------------------------------------------")
    print('==> TLA-setup finished. Total time elapsed: ', trun, '[HH:MM:SS]')
    plt.close('all')
    
    # %% end
    return (0)


# %% Argument parser
if __name__ == "__main__":

    # Create the parser
    my_parser = ArgumentParser(prog="tla_setup_sample",
                               description="# Single Sample Pre-processing " +
                               "module for Tumor Landscape Analysis #",
                               allow_abbrev=False)

    # Add the arguments
    my_parser.add_argument('-v', '--version',
                           action='version',
                           version="%(prog)s " + __version__)

    my_parser.add_argument('argsfile',
                           metavar="argsfile",
                           type=str,
                           help="Argument table file (.csv) for study set")

    my_parser.add_argument('casenum',
                           metavar="casenum",
                           type=int,
                           help="Set case number to be processed (zero based)")

    my_parser.add_argument("--graph",
                           default=False,
                           action="store_true",
                           help="If <<--graph>> is used, then print graphs")

    my_parser.add_argument("--redo",
                           default=False,
                           action="store_true",
                           help="If <<--redo>> is used, then redo analysis")

    # passes arguments object to main
    main(my_parser.parse_args())
