from astroquery.irsa import Irsa
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from ipywidgets import interact, fixed


def f_to_mag(flux):
    result = 22.5 - 2.5 * np.log10(flux)
    return result


def cmd_setup(ra_min, ra_max, dec_min, dec_max, GridSize=100, height=8, width=20):
   
    """unWISE sources scatter plot"""
    Irsa.ROW_LIMIT = 100000
    data_table = Irsa.query_region("m33", catalog="unwise_2019", spatial="Polygon",
            polygon=[SkyCoord(ra=ra_min,dec=dec_min,unit=(u.deg,u.deg),frame='icrs'),
                     SkyCoord(ra=ra_max,dec=dec_min,unit=(u.deg,u.deg),frame='icrs'),
                     SkyCoord(ra=ra_max,dec=dec_max,unit=(u.deg,u.deg),frame='icrs'),
                     SkyCoord(ra=ra_min,dec=dec_max,unit=(u.deg,u.deg),frame='icrs')])
    
    data_table['W1'] = f_to_mag(data_table['flux_1'])
    data_table['W2'] = f_to_mag(data_table['flux_2'])
    data_table['W2 - W1'] = data_table['W2'] - data_table['W1']
    
    fig, (ax1, ax2) = plt.subplots(1,2)
    fig.set_figheight(height)
    fig.set_figwidth(width)
    
    ax1.set_xlabel('W2 - W1', fontsize=20)
    ax1.set_ylabel('W2', fontsize=20)
    ax1.scatter(data_table['W2 - W1'],data_table['W2'], s=10,c='black')
    ax1.set_ylim(ax1.get_ylim()[::-1])
    ax1.set_title('unWISE', size=20)
    
    """HST sources 2d histogram"""
    hst_table = pd.read_csv('./small_data.csv')
    selection = hst_table.where((hst_table['RA'] > ra_min) & (hst_table['RA'] < ra_max) 
                                & (hst_table['DEC'] > dec_min) & (hst_table['DEC'] < dec_max))
    
    
    color = selection['F110W_VEGA'] - selection['F160W_VEGA']
    mag = selection['F160W_VEGA']
    ax2.hexbin(color, mag, gridsize=GridSize, bins='log')
    ax2.invert_yaxis()
    ax2.invert_xaxis()
    ax2.set_xlabel('F110W - F160W', fontsize=20)
    ax2.set_ylabel('F160W', fontsize=20)
    ax2.set_title('HST', fontsize=20)
    
    
def cmds(ra_min, ra_max, dec_min, dec_max):
    interact(cmd_setup, ra_min=fixed(ra_min), ra_max=fixed(ra_max),
             dec_min=fixed(dec_min), dec_max=fixed(dec_max), 
             GridSize=(50,600,1), height=(1,30,0.5), width=(1,30,0.5)
            )
