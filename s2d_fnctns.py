#!/usr/bin/env python

#import standard libraries
import obspy.imaging.beachball
import datetime
import os
import csv
import pandas as pd
import numpy as np
from geopy.distance import vincenty
from math import *
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


def stoc(input_file, output_file):

    '''Converts ssv and tsv files to csv'''
    
    with open(input_file) as fin, open(output_file, 'w') as fout:
        writer = csv.writer(fout)
        for line in fin:
            writer.writerow(line.split())


def castfloats(data):

    '''Casts all numerical and nan values to floats to avoid error in calculations'''
    
    data[['lat']] = data[['lat']].astype(float)
    data[['lon']] = data[['lon']].astype(float)
    data[['depth']] = data[['depth']].astype(float)
    data[['unc']] = data[['unc']].astype(float)
    if 'mag' in data.columns:
        data[['mag']] = data[['mag']].astype(float)
    if 'mrr' in data.columns:
        data[['mrr']] = data[['mrr']].astype(float)
        data[['mtt']] = data[['mtt']].astype(float)
        data[['mpp']] = data[['mpp']].astype(float)
        data[['mrt']] = data[['mrt']].astype(float)
        data[['mrp']] = data[['mrp']].astype(float)
        data[['mtp']] = data[['mtp']].astype(float)
    if 'Paz' in data.columns and 'Ppl' in data.columns:
        data[['Paz']] = data[['Paz']].astype(float)
        data[['Ppl']] = data[['Ppl']].astype(float)
        data[['Taz']] = data[['Taz']].astype(float)
        data[['Tpl']] = data[['Tpl']].astype(float)
        data[['S1']] = data[['S1']].astype(float)
        data[['D1']] = data[['D1']].astype(float)
        data[['R1']] = data[['R1']].astype(float)
        data[['S2']] = data[['S2']].astype(float)
        data[['D2']] = data[['D2']].astype(float)
        data[['R2']] = data[['R2']].astype(float)
    if 'strike' in data.columns and 'dip' in data.columns:
        data[['strike']] = data[['strike']].astype(float)
        data[['dip']] = data[['dip']].astype(float)
    return data


def rid_nans(df):

    '''Removes points where lat,lon,depth, or uncertainty values are not provided.'''
    
    df = df[np.isfinite(df['lat'])]
    df = df[np.isfinite(df['lon'])]
    df = df[np.isfinite(df['depth'])]
    df = df[np.isfinite(df['unc'])]
    return df


def write_data(df, output_file):

    ''' Arguments:  df - filtered dataframe to be written to file
                    output_file - output file where data is to be written  '''
    
    # If file name does not exist, creates file and writes filtered dataframe to it
    df = castfloats(df)
    df = rid_nans(df)
    if not os.path.isfile(output_file):
        with open(output_file, 'w') as f:
            df.to_csv(f, header=True, index=False, float_format='%0.3f', na_rep = float('nan'))

    # If the output file already exists, new filtered data points are appended to 
    # existing information
    else:
        old = pd.read_csv(output_file)
        all = pd.concat([old,df])
        all = castfloats(all)
        all = rid_nans(all)
        if len(df.columns) > len(old.columns):
            all = all[df.columns]
        else:
            all = all[old.columns]

        # Writes desired columns of a filtered dataframe to the output file
        with open(output_file, 'w') as f:
            all.to_csv(f, header=True, index=False, float_format='%0.3f', na_rep = float('nan'))


def inbounds(args, data):

    ''' Arguments:  args - input arguments provided from command line arguments
                    data - dataframe to be filtered based on bounds
                    
        Returns:    data - filtered dataframe based on bounds        '''
    
    # Eliminates data points that are not within specified bounds where provided
    if 'time' in data.columns:
        data['time'] = pd.to_datetime(data['time'])
    
    stime = datetime.datetime(1900,1,1)
    etime = datetime.datetime.utcnow()
    
    if args.startTime and args.endTime and args.startTime >= args.endTime:
        print 'End time must be greater than start time.  Your inputs: Start %s \
        End %s' % (args.startTime, args.endTime)
        sys.exit(1)

    if args.bounds is not None:
        lonmin = args.bounds[0]
        lonmax = args.bounds[1]
        latmin = args.bounds[2]
        latmax = args.bounds[3]
        minwest = lonmin > 0 and lonmin < 180
        maxeast = lonmax < 0 and lonmax > -180
        if minwest and maxeast:
            data = data[(data.lon >= lonmin) | (data.lon <= lonmax)]
        else:
            data = data[(data.lon >= lonmin) & (data.lon <= lonmax)]
        data = data[(data.lat >= latmin) & (data.lat <= latmax)]

    if args.startTime is not None and 'time' in data.columns:
        stime = args.startTime
        data = data[data.time >= stime]

    if args.endTime is not None and 'time' in data.columns:
        etime = args.endTime
        data = data[data.time <= etime]

    if args.magRange is not None and 'mag' in data.columns:
        magmin = args.magRange[0]
        magmax = args.magRange[1]
        data = data[(data.mag >= magmin) & (data.mag <= magmax)]

    return data


def cmtfilter(data):

    ''' Arguments:  data - data with all shallow/nonshallow and thrust/nonthrust earthquake
    
        Returns:    filtered - fitered dataframe which DEPENDS ON WHAT YOU DO/DONT COMMENT OUT
        
                                (1) filters only shallow earthquakes that have MT criteria which are non thrust
                                all other shallow earthquakes WITHOUT MT info are NOT filtered
                                
                                OR
                                
                                (2) filters ALL shallow earthquakes UNLESS they have MT info and that
                                MT info has the criteria of a thrust event. '''
    
    # Removes non-thrust events from depths shallower than 60 km
    deep_data = data[data.depth >= 70]
    
    # Includes shallow data without MT info (1) - comment out next two lines for (2)
    dfn = data[np.isnan(data['Paz'])]
    dfn = dfn[data.depth < 70]
    
    data = data[np.isfinite(data['Paz'])]
    shallow_data = data[data.depth < 70]
    
    # Depending on which MT info are provided, filters non-thrust, shallow events
    if 'Ndip' in shallow_data.columns:
        thrust_rake = (shallow_data.Tpl>50) & (shallow_data.Ndip<=30)
    else:
        thrust_rake = ((shallow_data.R1>30) & (shallow_data.R2>30)
                       & (shallow_data.R1<150) & (shallow_data.R2<150))
    shallow_data = shallow_data[thrust_rake]
    
    # Includes shallow data without MT info (1) - comment out next line for (2)
    filtered = pd.concat([deep_data, shallow_data, dfn])

    # Only includes shallow thrust events (2) - uncomment line below for (2) and comment necessary lines above
    # filtered = pd.concat([deep_data, shallow_data])
    
    # Rearranges columns / filters out unecessary columns
    filtered=filtered[['lat','lon','depth','unc','etype','mag','time',
                       'Paz','Ppl','Taz','Tpl','S1','D1','R1','S2','D2','R2',
                       'strike','dip']]
    return filtered


def make_moment_tensor(mrr,mtt,mpp,mrt,mrp,mtp): #r,t,p = x,y,z

    '''Used in m_to_planes below. Makes a moment tensor object from moment tensor components'''
    
    return obspy.imaging.beachball.MomentTensor(mrr,mtt,mpp,mrt,mrp,mtp,1)


def m_to_planes(mrr,mtt,mpp,mrt,mrp,mtp,n):

    '''Takes a moment tensor and calculates the P, N, and T axes and nodal plane information.
        
        Used in moment_calc below.  Returns one of these values as specified by input (n). 
        The integer input specifies which index of the array of outputs to return.  '''
    
    mt = make_moment_tensor(mrr,mtt,mpp,mrt,mrp,mtp)
    #axes = obspy.imaging.beachball.MT2Axes(mt) #returns T, N, P
    #fplane = obspy.imaging.beachball.MT2Plane(mt)#returns strike, dip, rake
    #aplane = obspy.imaging.beachball.AuxPlane(fplane.strike, fplane.dip, fplane.rake)
    #MAF changed because functions use lowercase, and aux_plane name includes underscore
    axes = obspy.imaging.beachball.mt2axes(mt) #returns T, N, P
    fplane = obspy.imaging.beachball.mt2plane(mt)#returns strike, dip, rake
    aplane = obspy.imaging.beachball.aux_plane(fplane.strike, fplane.dip, fplane.rake)
    Tstrike = axes[0].strike
    Tdip = axes[0].dip
    Nstrike = axes[1].strike
    Ndip = axes[1].dip
    Pstrike = axes[2].strike
    Pdip = axes[2].dip
    S1 = fplane.strike
    D1 = fplane.dip
    R1 = fplane.rake
    S2 = aplane[0]
    D2 = aplane[1]
    R2 = aplane[2]
    mplanes = [Pstrike,Pdip,Tstrike,Tdip,S1,D1,R1,S2,D2,R2,Nstrike,Ndip]
    return mplanes[n]

def moment_calc(df, args):

    ''' Creates and appends columns with Principal Axis and Nodal Plane information.
        
        Used in makeframe below. Takes moment tensor information from input dataframe 
        columns and creates 11 new columns with information used to distinguish between thrust
        and non-thrust earthquakes.   
        
        Arguments:  df - dataframe with mt information in the form mrr,mtt,mpp,mrt,mrp,mtp
                    args - input arguments provided from command line arguments   
                    
        Returns:    df - dataframe with mt information in the form Paz,Ppl,Taz,Tpl,S1,D1,R1,S2,D2,R2 '''
    
    
    try:
        # Only calculates MT info where it exists in EQ datasets
        df = inbounds(args, df)
        dfm = df[np.isfinite(df['mrr'])]
        dfn = df[df['mrr'].isnull()]
    except:
        raise Exception,'If file contains earthquake information (event-type = EQ), \
        required columns include: lat,lon,depth,mag,time. The columns of the current \
        file: %s. Check file format to ensure these columns are present and properly \
        labeled.' % df.columns

    # Calculates each new column of MT info
    try:
        dfm['Paz']=dfm.apply(lambda row: m_to_planes(row['mrr'], row['mtt'], row['mpp'],
                                                     row['mrt'], row['mrp'], row['mtp'],0),axis=1)
        dfm['Ppl']=dfm.apply(lambda row: m_to_planes(row['mrr'], row['mtt'], row['mpp'],
                                                     row['mrt'], row['mrp'], row['mtp'],1),axis=1)
        dfm['Taz']=dfm.apply(lambda row: m_to_planes(row['mrr'], row['mtt'], row['mpp'],
                                                     row['mrt'], row['mrp'], row['mtp'],2),axis=1)
        dfm['Tpl']=dfm.apply(lambda row: m_to_planes(row['mrr'], row['mtt'], row['mpp'],
                                                     row['mrt'], row['mrp'], row['mtp'],3),axis=1)
        dfm['S1']=dfm.apply(lambda row: m_to_planes(row['mrr'], row['mtt'], row['mpp'],
                                                    row['mrt'], row['mrp'], row['mtp'],4),axis=1)
        dfm['D1']=dfm.apply(lambda row: m_to_planes(row['mrr'], row['mtt'], row['mpp'],
                                                    row['mrt'], row['mrp'], row['mtp'],5),axis=1)
        dfm['R1']=dfm.apply(lambda row: m_to_planes(row['mrr'], row['mtt'], row['mpp'],
                                                    row['mrt'], row['mrp'], row['mtp'],6),axis=1)
        dfm['S2']=dfm.apply(lambda row: m_to_planes(row['mrr'], row['mtt'], row['mpp'],
                                                    row['mrt'], row['mrp'], row['mtp'],7),axis=1)
        dfm['D2']=dfm.apply(lambda row: m_to_planes(row['mrr'], row['mtt'], row['mpp'],
                                                    row['mrt'], row['mrp'], row['mtp'],8),axis=1)
        dfm['R2']=dfm.apply(lambda row: m_to_planes(row['mrr'], row['mtt'], row['mpp'],
                                                    row['mrt'], row['mrp'], row['mtp'],9),axis=1)
        dfm['Ndip']=dfm.apply(lambda row: m_to_planes(row['mrr'], row['mtt'], row['mpp'],
                                                      row['mrt'], row['mrp'], row['mtp'],11),axis=1)
                                        
        # Concatenates events with and without MT info
        dfm = cmtfilter(dfm)
        df = pd.concat([dfm,dfn])
        
        # Rearranges columns and returns
        df = df[['lat','lon','depth','unc','etype','mag','time',
               'Paz','Ppl','Taz','Tpl','S1','D1','R1','S2','D2','R2',
               'strike','dip']]
        return df
    except:
    
        # if exception is caught, try to return only events without MT info
       try:
           if len(dfm) == 0:
               return dfn
       except:
           raise Exception,'Where moment tensor information is available, columns \
           must be labeled: mrr,mpp,mtt,mrp,mrt,mtp'

def ymdhmsparse(input_file):

    '''Parses Yr Mo Day Hr Min Sec into one datetime object when provided in distinguished columns.
        
        Used in makeframe below. Returns a new dataframe with parsed datetimes.    '''
    
    ymdhms = {'time':['year','month','day','hour','min','sec']}
    dparse = lambda x: pd.datetime.strptime(x, '%Y %m %d %H %M %S')
    cols = ['year','month','day','hour','min','sec','lat','lon','depth','mag']
    data = pd.read_csv(input_file, parse_dates=ymdhms, usecols=cols, date_parser=dparse)
    return data

def raiseUnc(x):

    ''' Raises unreasonably low uncertainties for earthquakes to a value greater 
        than that of average active source data points (which is 5 km).   '''
    
    if x < 6:
        return 6
    else:
        return x

def makeframe(data, fcsv, event_type, uncertainty, args):

    ''' Arguments:  data - semi-filtered data frame to be filtered more and written to file
                    fcsv - filename of output file
                    event_type - kind of data i.e. BA, EQ, ER, TO etc
                    uncertainty - unc value provided in command line or set by default for etype
                    args - input arguments provided from command line arguments 
                    
        Returns:    data - fully filtered dataset to be written to output file '''
    
    # Parses Yr Mo Day Hr Min Sec into one datetime object when provided in distinguished columns
    if 'year' in data.columns and 'sec' in data.columns and 'mag' in data.columns:
        data = ymdhmsparse(fcsv)
    
    # If ISC-GEM data is provided, high quality, low uncertainties are included in place of
    # the default values assigned in s2d.py main method.
    if 'unc' in data.columns and 'q' in data.columns:
        try:
            data = data[(data.uq != 'C') & (data.unc < uncertainty)]
        except:
            raise Exception,'When adding a file with uncertainty quality, the column \
            representing that quality must be labeled as uq'

    # uses OG uncertainties where provided. Raises them if they are unreasonably low
    elif 'unc' in data.columns:
        uncert = data['unc'].values
        if isnan(uncert[1]):
            data['unc'] = uncertainty
        elif event_type == 'EQ':
            data['unc'] = data.apply(lambda row: raiseUnc(row['unc']),axis=1)
        else:
            pass

    # If no uncertainty column is included, the one provided in command line arguments is
    # used to add a new column to the data, alternatively, the default value assigned in s2d.py is used
    else:
        data['unc'] = uncertainty
    pd.options.mode.chained_assignment = None

    # A new column marking the event type is added to the data. Everything is cast as a float
    data['etype'] = event_type
    data = castfloats(data)

    # Adds strike and dip from tomography data where available
    if 'strike' in data.columns and 'dip' in data.columns and event_type == 'TO':
        data = data[['lat','lon','depth','unc','etype','strike','dip']]
    else:
        data['strike'] = float('nan')
        data['dip'] = float('nan')
    if len(data.columns) == 7:
        data = data[['lat','lon','depth','unc','etype','strike','dip']]

    # Calculates moment tensor info where applicable and removes shallow, non-thrust events
    if 'mrr' in data.columns:
        data = moment_calc(data,  args)
    elif 'Paz' in data.columns and event_type == 'EQ':
        data = data[['lat','lon','depth','unc','etype','mag','time',
                     'Paz','Ppl','Taz','Tpl','S1','D1','R1','S2','D2','R2',
                     'strike','dip']]
        data = cmtfilter(data)
    elif 'time' in data.columns and 'mag' in data.columns:
        data = data[['lat','lon','depth','unc','etype','mag','time','strike','dip']]
    else:
        pass

    return data

##########################################################################################################
#The following serves to create a rough plot of the data types compiled with s2d.py.
##########################################################################################################

def plot_map(lons, lats, c, legend_label, projection='mill',
             llcrnrlat=-80, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180, resolution='i'):
    
    ''' Optional Arguments: projection - map projection, default set as 'mill'
                            llcrnrlat - lower left corner latitude value, default is -80
                            urcrnrlat - upper right corner latitude value, default is 90
                            llcrnrlon - lower left corner longitude value, default is -180
                            urcrnrlon - upper right corner longitude value, default is 180
                            resolution - the resolution of the plot, default is 'i'
                            
        Required Arguments: lons - list of longitude values to be plotted
                            lats - list of latitude values to be plotted
                            c - the color of the points to be plotted
                            legend_label - how this set of points will be labeled on the legend
                            
        Returns:            m - a basemap object defined by input bounds with input points included '''
    
    # Creates a basic plot of a series of lat,lon points over a defined region
    m = Basemap(projection=projection, llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,
                llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon, resolution=resolution)
    m.drawcoastlines()
    m.drawmapboundary()
    m.drawcountries()
    m.etopo()
    m.drawmeridians(np.arange(llcrnrlon, urcrnrlon, 5), labels=[0,0,0,1], fontsize=10)
    m.drawparallels(np.arange(llcrnrlat, urcrnrlat, 5), labels=[1,0,0,0], fontsize=10)
    x,y = m(lons, lats)
    m.scatter(x, y, color=c, label=legend_label, marker='o', edgecolor='none', s=10)
    return m

def datelinecross(x):

    '''Converts negative longitudes to their positive equivalent for the sake of plotting.'''
    
    if x<0:
        return x+360
    else:
        return x

def slabplotter(args):

    '''
        '''
    
    # Gathers events and bounds from input arguments
    file = args.outFile
    lonmin = args.bounds[0]
    lonmax = args.bounds[1]
    latmin = args.bounds[2]
    latmax = args.bounds[3]
    data = pd.read_csv(file, low_memory=False)
    data = castfloats(data)
    data = data[['lat','lon','depth','etype','Paz']]
    xdim = (abs(lonmax-lonmin)) * .01
    ydim = (abs(latmax-latmin)) * .01

    # Handles input files that cross the dateline
    minwest = lonmin > 0 and lonmin < 180
    maxeast = lonmax < 0 and lonmax > -180
    if minwest and maxeast:
        data['lon'] = data.apply(lambda row: datelinecross(row['lon']),axis=1)
        lonmax += 360
        xdim = ((180-lonmin) + (180+lonmax)) *.01

    # Distinguishes each data point by event type
    dataCMT = data[np.isfinite(data['Paz'])]
    data = data[pd.isnull(data).any(axis=1)]
    dataEQ = data[data.etype == 'EQ']
    dataER = data[data.etype == 'ER']
    dataRF = data[data.etype == 'RF']
    dataTO = data[data.etype == 'TO']
    dataBA = data[data.etype == 'BA']
    dataAS = data[data.etype == 'AS']
    dataMT = data[data.etype == 'MT']

    # Creates a tuple for each event type with a list of lons, lats, and legend label
    dataEQt = (dataEQ['lon'].values, dataEQ['lat'].values,
               'EQ: %i points' % len(dataEQ))
    dataCMTt = (dataCMT['lon'].values, dataCMT['lat'].values,
                'EQ with MT: %i points' % len(dataCMT))
    dataERt = (dataER['lon'].values, dataER['lat'].values,
               'ER: %i points' % len(dataER))
    dataRFt = (dataRF['lon'].values, dataRF['lat'].values,
               'RF: %i points' % len(dataRF))
    dataTOt = (dataTO['lon'].values, dataTO['lat'].values,
               'TO: %i points' % len(dataTO))
    dataASt = (dataAS['lon'].values, dataAS['lat'].values,
               'AS: %i points' % len(dataAS))
    dataBAt = (dataBA['lon'].values, dataBA['lat'].values,
               'BA: %i points' % len(dataBA))
    dataMTt = (dataMT['lon'].values, dataMT['lat'].values,
               'MT: %i points' % len(dataMT))

    # Makes a list of tuples made above and associated colors to be plotted
    typelist = [dataBAt,dataEQt,dataCMTt,dataERt,dataRFt,dataTOt,dataASt]
    colorlist = ['black','white','gold','red','blueviolet','turquoise',
                 'lime','darkorange','blue','darkcyan','orchid']

    # Makes the image an appropriate size relative to the bound dimensions
    while xdim < 14:
        xdim *= 1.1
    while ydim < 8:
        ydim *= 1.1
    fig = plt.figure(figsize=(xdim,ydim))
    ax = plt.subplot(111)
    color = 0

    # Plots the amount and type of data in a given file over a designated region.
    for type in typelist:
        if len(type[0]) > 0:
            plot_map(type[0], type[1], colorlist[color], type[2],
                     llcrnrlat=latmin, urcrnrlat=latmax, llcrnrlon=lonmin, urcrnrlon=lonmax)
            color += 1
    ax.set_title(file[:-4])
    lgd = ax.legend(loc='center left', bbox_to_anchor=(1,0.5), fontsize=12)

    # The plot is saved to the current directory
    try:
        fig.savefig(file[:-4]+'types.png', format='png', bbox_extra_artists=(lgd,),
                    bbox_inches='tight', pad_inches=0.5)
    except:
        raise Exception,'No information was written to the file, check bounds to ensure that \
        they encompass subduction zone information within the bounds of the input file(s)'

    plt.show()


##############################################################################################
#Everything below this point serves the purpose of identifying and
#eliminating duplicate events between multiple earthquake catalog entries.
##############################################################################################

class Earthquake:

    '''Creates an earthquake object from which event information can be extracted'''
    
    def __init__(self,time,coords,depth,lat,lon,mag,catalog):
        self.time = time
        self.coords = coords
        self.depth = depth
        self.lat = lat
        self.lon = lon
        self.mag = mag
        self.catalog = catalog

def getvals(row):

    '''Gathers time, lat, lon, depth, mag, information from row in dataframe.'''
    
    time = row['time']
    lat = row['lat']
    lon = row['lon']
    depth = row['depth']
    mag = row['mag']
    ep = (lat,lon)
    return time,ep,depth,lat,lon,mag

def boundtrim(cat1, cat2):

    ''' Arguments:  cat1 - an earthquake catalog to be compared with cat2
                    cat2 - an earthquake catalog to be compared to cat1
                    
        Returns:     cat1, cat2 - trimmed earthquake catalogs that only extend across bounds
                                    where they both exist. Reduces processing time
        '''
    
    # Trims two earthquake catalogs to fit over the same region
    lonmin1, lonmin2 = cat1['lon'].min(), cat2['lon'].min()
    latmin1, latmin2 = cat1['lat'].min(), cat2['lat'].min()
    lonmax1, lonmax2 = cat1['lon'].max(), cat2['lon'].max()
    latmax1, latmax2 = cat1['lat'].max(), cat2['lat'].max()
    minwest = (lonmax1 > 0 and lonmax1 < 180) or (lonmax2 > 0 and lonmax2 < 180)
    maxeast = (lonmin1 < 0 and lonmin1 > -180) or (lonmin2 < 0 and lonmin2 > -180)
    difference = abs(lonmin1-lonmax1)>180 or abs(lonmin2-lonmax2)>180
    if minwest and maxeast and difference:
        pass
    else:
        cat1 = cat1[(cat1.lon >= lonmin2) & (cat1.lon <= lonmax2)]
        cat2 = cat2[(cat2.lon >= lonmin1) & (cat2.lon <= lonmax1)]
    cat1 = cat1[(cat1.lat >= latmin2) & (cat1.lat <= latmax2)]
    cat2 = cat2[(cat2.lat >= latmin1) & (cat2.lat <= latmax1)]
    return cat1, cat2

def timetrim(cat1, cat2):

    ''' Arguments:  cat1 - an earthquake catalog to be compared with cat2
                    cat2 - an earthquake catalog to be compared to cat1
                    
        Returns:     cat1, cat2 - trimmed earthquake catalogs that only extend across time
                                    frames where they both exist. Reduces processing time
        '''
    
    # Trims two earthquake catalogs to fit over the same time range
    cat1['time'] = pd.to_datetime(cat1['time'])
    cat2['time'] = pd.to_datetime(cat2['time'])
    cat1min, cat1max = cat1['time'].min(), cat1['time'].max()
    cat2min, cat2max = cat2['time'].min(), cat2['time'].max()
    cat1 = cat1[(cat1.time >= cat2min) & (cat1.time <= cat2max)]
    cat2 = cat2[(cat2.time >= cat1min) & (cat2.time <= cat1max)]
    return cat1, cat2

def earthquake_string(eqo):

    ''' Puts earthquake information into a string to be written or printed
    
        Arguments:  eqo - earthquake object
        
        Returns:    eqos - a string of information stored in earthquake object input argument '''
    
    eqos = (str(eqo.lat) + ',' + str(eqo.lon) + ',' + str(eqo.depth) + ','
            + str(eqo.mag) + ',' + str(eqo.time) + ',' + eqo.catalog)
    return eqos

def find_closest(eqo, eqm1, eqm2):

    '''Determines which of two potential matches in one catalog is closer to an event in another.
        
        Arguments:  eqo - earthquake event in first catalog that matches two events in the second
                    eqm1 - the first event in the second catalog that matches eqo
                    eqm2 - the second event in the second catalog that matches eqo
                    
        Returns:    closest - the closest event weighting time first, then distance, then magnitude '''
    
    # Prints information to console to make user aware of more than one match
    print '-------------------------------------- lat %s lon %s depth %s mag %s time \
        %s catlog' % (',',',',',',',',',')
    print 'There is more than one match for event: %s' % earthquake_string(eqo)
    print 'event1: %s' % earthquake_string(eqm1)
    print 'event2: %s' % earthquake_string(eqm2)
    
    # Gets distance between either event and the common match eqo
    darc1 = vincenty(eqo.coords, eqm1.coords).meters/1000
    darc2 = vincenty(eqo.coords, eqm2.coords).meters/1000
    dh1 = abs(eqo.depth - eqm1.depth)
    dh2 = abs(eqo.depth - eqm2.depth)
    dist1 = sqrt(darc1*darc1 + dh1*dh1)
    dist2 = sqrt(darc2*darc2 + dh2*dh2)
    
    # Gets magnitude and time differences between each event and the common match
    dtime1 = abs(eqo.time - eqm1.time)
    dtime2 = abs(eqo.time - eqm2.time)
    dmag1 = abs(eqo.mag - eqm1.mag)
    dmag2 = abs(eqo.mag - eqm2.mag)
    
    # Finds the closest match to eqo by checking time first, then distance, then magnitude
    if dtime1 < dtime2:
        closest = eqm1
    elif dtime2 < dtime1:
        closest = eqm2
    elif dtime1 == dtime2 and dist1 < dist2:
        closest = eqm1
    elif dtime1 == dtime2 and dist2 < dist1:
        closest = eqm1
    elif dmag1 == dmag2 and dist1 == dist2 and dmag1 < dmag2:
        closest = eqm1
    elif dmag1 == dmag2 and dist1 == dist2 and dmag2 < dmag1:
        closest  = eqm2

    # If all things are equal, the first event is chosen as a match by default
    else:
        print 'The two events are equidistant to the match in time, space, and magnitude.\
            The second event was therefore determined independent.'
        closest = eqm1
        return closest
    print '>>>>closest event: %s' % earthquake_string(closest)
    return closest

def removematches(dfo, dfm):

    '''Eliminates events in dfo (dataframe) that are found in dfm (dataframe) '''
    
    ind = (dfo.time.isin(dfm.time) & dfo.lat.isin(dfm.lat) & dfo.lon.isin(dfm.lon)
           & dfo.mag.isin(dfm.mag) & dfo.depth.isin(dfm.depth))
    dfo = dfo[~ind]
    return dfo

def rid_matches(cat1, cat2, name1, name2):

    ''' Compares two catalogs, identifies and removes matching events from cat2.
        
        Arguments:  cat1 - the first catalog (dataframe), no events are removed from this catalog
                    cat2 - the second catalog (dataframe), events in this catalog that are close 
                            in space, time, and magnitude to those in cat1 are filtered out  
                    name1 - the name of the first catalog, used for printing/bookeeping purposes
                    name2 - the name of the second catalog, used for printing/bookeeping purposes
                    
        Returns:    df - a filtered version of cat2 without events that match those in cat1 '''
    
    # Setting constants that define matching criteria
    tdelta = 30
    distdelta = 100
    magdelta = 0.5
    
    # Ensuring that all times are in datetime object format & trimming catalogs to only extend
    # accross the bounds and time constraints of the other
    cat1['time'] = pd.to_datetime(cat1['time'])
    cat2['time'] = pd.to_datetime(cat2['time'])
    cat1c,cat2c = timetrim(cat1, cat2)
    cat1c,cat2c = boundtrim(cat1c, cat2c)
    
    # Making dataframe/filename to store matching events for bookeeping
    try:
        name1w = name1[-10:]
        name2w = name2[-10:]
    except:
        name1w = name1[:-4]
        name2w = name2[:-4]
    matches = pd.DataFrame(columns = ['lat','lon','depth','mag','time','catalog'])
    count = 0

    # Compares each event in cat2 to each event in cat1
    for index,row in cat1c.iterrows():
        n = 0
        
        # Getting earthquake info from event and storing it in an Earthquake object (cat1)
        time1, ep1, depth1, lat1, lon1, mag1 = getvals(row)
        eq1 = Earthquake(time1, ep1, depth1, lat1, lon1, mag1, name1w)
        for index, r in cat2c.iterrows():
        
            # Getting earthquake info from event and storing it in an Earthquake object (cat1)
            time2, ep2, depth2, lat2, lon2, mag2 = getvals(r)
            eq2 = Earthquake(time2, ep2, depth2, lat2, lon2, mag2, name2w)
            
            # If events are close in time, space, and magnitude add event from cat2 to match list
            if abs(time1-time2) < datetime.timedelta(seconds = tdelta):
                if vincenty(ep1,ep2).meters/1000 <= distdelta:
                    if abs(mag1-mag2) < magdelta:
                        
                        # If there is already a match for this event, find the closest
                        # The closest is stored and compared to third, fourth matches etc if they exist
                        if n >= 1:
                            match = find_closest(eq1, match, eq2)
                            n += 1
                        if n == 0:
                            match = eq2
                            n += 1
                    else:
                        pass
                else:
                    pass
            else:
                pass
        
        # Add matching events to match dataframe
        if n > 0:
            lat1,lon1,depth1,mag1,time1,name1
            matches.loc[len(matches)+1] = [lat1, lon1, depth1, mag1, time1, name1w]
            matches.loc[len(matches)+1] = [match.lat, match.lon, match.depth, match.mag,
                                           match.time, name2w]
            count += 1
        
    # Write matches to matching file
    matchfile = name1w + name2w + '-matches.csv'
    
    # Remove matches from cat2
    df = removematches(cat2, matches)
    
    # Print general results to console
    print '%i matches were found between the catalogs: %s and %s.' % (count, name1, name2)
    if count > 0:
        with open(matchfile,'w') as f:
            matches.to_csv(f, header=True, index=False, float_format='%0.4f')
        print 'The pairs can be found in the file: ** %s **, which has been written added to the current directory.' % (name1w + name2w + '-matches.csv')
        print 'Based on the order of entry, in the instance of duplicate events, the entries in ** %s ** were added to the slab file while the entries in ** %s ** were not added.' % (name1, name2)

    # Return filtered catalog to be written to output file
    return df
