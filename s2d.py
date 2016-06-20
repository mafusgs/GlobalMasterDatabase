#!/usr/bin/env python

#stdlib
import argparse
from datetime import datetime,timedelta
import os
import sys
import re
import pandas as pd
import numpy as np

#local library
import s2d_fnctns

class NewFile:

    '''Creates a file object with associated uncertainty and event type'''
    
    def __init__(self, filename, unc, event_type):
        self.filename = filename
        self.event_type = event_type
        self.unc = unc

def maketime(timestring):

    '''Used in argument parser below. Makes a datetime object from a timestring.'''
    
    TIMEFMT = '%Y-%m-%dT%H:%M:%S' 
    #TIMEFMT = ' %Y-%m-%d %H:%M:%S' #does this need to change?
    DATEFMT = '%Y-%m-%d'
    outtime = None
    try:
        outtime = datetime.strptime(timestring, TIMEFMT)
    except:
        try:
            outtime = datetime.strptime(timestring, DATEFMT)
        except:
            raise Exception,'Could not parse time or date from %s' % timestring
    return outtime

def infile(s):

    '''Stores filename, event type, and uncertainty where provided from comma separated string.'''
    
    default_uncertainty = 15
    try:
        infile,unc,etype = s.split(',')
        unc = float(unc)
        return (infile, unc, etype)
    except:
        try:
            s = s.split(',')
            infile, unc, etype = s[0], default_uncertainty, s[1]
            return (infile, unc, etype)
        except:
            raise argparse.ArgumentTypeError('Input file information must be \
                                             given as infile,unc,etype or as infile,etype')

def writetofile(input_file, output_file, event_type, uncertainty, args, catalogs, file_no):

    ''' Writes an input file object to the given output file.
        
        Acquires the necessary columns from the file, calculates moment tensor information.
        Eliminates rows of data that do not fall within the specified bounds
            (date, magnitude, & location).
        If the event type is an earthquake, the catalog is compared to all previously
            entered catalogs. Duplicate events are removed from the subsequent entries
            (prioritization is determined by the order in which catalogs are entered).
        Writes filtered dataframe to output file and prints progress to console.
        
        Arguments:  input_file - input file from input or slab2database
                    output_file - file where new dataset will be written
                    event_type - two letter ID that indicates the type of data (AS, EQ, BA, etc)
                    uncertainty - the default uncertainty associated with this file or event type
                    args - arguments provided from command line (bounds, magnitude limits, etc)
                    catalogs - a list of EQ catalogs that are being written to this file
                    file_no - file number, used for making event IDs '''
    
    in_file = open(input_file)
    fcsv = (input_file[:-4]+'.csv')
    # Reading .csv file into dataframe - all files must be in .csv format
    try:
        if input_file.endswith('.csv'):
            data = pd.read_csv(input_file, low_memory=False)
        else:
            print 'Input file %s was not written to file. MUST BE IN .CSV FORMAT' % input_file
            pass
    except:
        raise Exception,'Could not read file %s. A header line of column labels \
        followed by a deliminated dataset is expected. Check file format to ensure this \
        is such. All files must be in .csv format.' % input_file
    data = s2d_fnctns.makeframe(data, fcsv, event_type, uncertainty, args)
    data = s2d_fnctns.inbounds(args, data)
    # Removing fixed events
    if event_type == 'EQ':
        data = data[(data.depth > 33.001) | (data.depth < 33.0)]
        data = data[(data.depth > 10.001) | (data.depth < 10.0)]
        data = data[(data.depth > 15.001) | (data.depth < 15.0)]
        data = data[(data.depth > 35.001) | (data.depth < 35.0)]
        data = data[(data.depth > 100.001) | (data.depth < 100.0)]
        data = data[(data.depth > 20.001) | (data.depth < 20.0)]
        data = data[(data.depth > 25.001) | (data.depth < 25)]
        data = data[(data.depth > 50.001) | (data.depth < 50.0)]
        data = data[(data.depth > 47.001) | (data.depth < 47.0)]
        data = data[(data.depth > 150.001) | (data.depth < 150.0)]
    #Removing duplicate entries for the same event
        try:
           tup = (data, fcsv)
           if len(catalogs) > 0: #i.e. if there are multiple catalogues...
                for cat in catalogs:
                    data = s2d_fnctns.rid_matches(cat[0], data, cat[1], fcsv) 
           elif len(catalogs) == 0:
                catalogs.append(tup)
        except:
            raise Exception,'If file contains earthquake information (event-type = EQ), \
            required columns include: lat,lon,depth,mag,time. The columns of the current \
            file: %s. Check file format to ensure these columns are present and properly \
            labeled.' % data.columns
    # Removing fixed events
    if event_type == 'ER':
        data = data[(data.depth > 33.001) | (data.depth < 33.0)]
        data = data[(data.depth > 10.001) | (data.depth < 10.0)]
        data = data[(data.depth > 15.001) | (data.depth < 15.0)]
        data = data[(data.depth > 35.001) | (data.depth < 35.0)]
        data = data[(data.depth > 100.001) | (data.depth < 100.0)]
        data = data[(data.depth > 47.001) | (data.depth < 47.0)]
        data = data[(data.depth > 150.001) | (data.depth < 150.0)]
    start_ID = file_no*100000
    stop_ID = start_ID + len(data)
    ID = np.arange(start_ID, stop_ID, 1)
    data['ID'] = ID
    s2d_fnctns.write_data(data, output_file)
    print 'The file: %s was written to %s' % (input_file, output_file)
    print '---------------------------------------------------------------------------------'

def getslab(args):

    ''' Arguments:  args - input arguments from command line
        
        Returns: slablist - list of slab regions that fall within the specified bounds  '''
    
    # Defines slab regions
    #(lonmin,lonmax,latmin,latmax,ID)
    alu = (165,-140,50,65,'alu')
    mex = (-106,-80,6,20,'mex')
    cas = (-128,-120,38,52,'cas')
    izu = (135,150,10,36,'izu')
    ker = (175,-170,-38,-14,'ker')
    kur = (129,166,31,57,'kur')
    phi = (123,128,7,15,'phi')
    ryu = (121,139,22,39,'ryu')
    van = (164,173,-24,-10,'van')
    sco = (-30,-22,-62,-55,'sco')
    sol = (145,163,-12,-3,'sol')
    sam = (-90,-60,-50,10,'sam')
    sum = (90,125,-12,12,'sum')
    slablist = [alu,mex,cas,izu,ker,kur,phi,ryu,van,sco,sol,sam,sum]
    
    # Compares bound arguments to specified slab ranges
    slabrange = []
    lonmin = args.bounds[0]
    lonmax = args.bounds[1]
    latmin = args.bounds[2]
    latmax = args.bounds[3]
    minwest = lonmin > 0 and lonmin < 180
    maxeast = lonmax < 0 and lonmax > -180
    
    # If the bounds cross the dateline
    if minwest and maxeast:
        for slab in slablist:
            lonrange = ((lonmax >= slab[0] or lonmax <= slab[1]) or
                        (lonmin >= slab[0] or lonmin <= slab[1]))
            latrange = ((latmin >= slab[2] and latmin <= slab[3]) or
                        (latmax <= slab[3] and latmax >= slab[2]))
            if lonrange and latrange:
                slabrange.append(slab[4])
    
    # If the bounds do not cross the dateline
    else:
        for slab in slablist:
            lonrange = ((lonmin >= slab[0] and lonmin <= slab[1]) or
                        (lonmax <= slab[1] and lonmax >= slab[0]))
            latrange = ((latmin >= slab[2] and latmin <= slab[3]) or
                        (latmax <= slab[3] and latmax >= slab[2]))
            if lonrange and latrange:
                slabrange.append(slab[4])
    
    # Return a list of slabs that fall within input bounds
    if len(slabrange) > 0:
        return slabrange
    else:
        print 'slab2database has no files for the specified region'

def main(args):

    '''What is executed upon running the program. 
        
        Runs through the provided input files and calls writetofile (above) for each. 
        The result is a new comma delimited file with information pertinent to Slab 2.0.  '''
    
    slab_ = args.outFile
    default_uncertainty = 15 #km
    
    filelist = []
    
    # If the database is called, compiles all files within bounds into a list to be written to
    # output file.
    if args.database is not None:
        for filename in os.listdir(args.database):
            if filename.endswith('.csv'):
                slabname,etype,name = filename.split('_')
                # Changes the default uncertainty based on event type.
                if etype == 'AS':
                    default_uncertainty = 2.5
                elif etype == 'RF':
                    default_uncertainty = 10
                elif etype == 'ER':
                    default_uncertainty = 10
                elif etype == 'BA':
                    default_uncertainty = 1
                # If bounds are provided, only adds files labeled with the associated slab region to file list.
                if args.bounds:
                    slablist = getslab(args)
                    if slabname in slablist:
                        # Creates file object to be written to output file
                        f = NewFile(args.database+'/'+filename, default_uncertainty, etype)
                        filelist.append(f)
                    elif slabname == 'ALL':
                        # Creates file object to be written to output file
                        f = NewFile(args.database+'/'+filename, default_uncertainty, etype)
                        filelist.append(f)
                    else:
                        pass
                else:
                    # Creates file object to be written to output file
                    f = NewFile(args.database+'/'+filename, default_uncertainty, etype)
                    filelist.append(f)
            else:
                print 'The file %s was not written to the dataset - FILES MUST BE IN .CSV FORMAT' % filename
    # Adds additional input files to file list if they are provided in command line.
    if args.input is not None:
        for file in args.input:
            try:
                f=NewFile(file[0], file[1], file[2])
                filelist.append(f)
            except:
                try:
                    f = NewFile(file[0], default_uncertainty, file[1])
                    filelist.append(f)
                except:
                    raise Exception,'Input file information must be given as infile,unc,etype \
                        or as infile,etype'
    # Notifies user and exits if no files are written to the output file.
    if len(filelist) == 0:
        print 'No input files were provided. Include argument: -i infile1,unc1,etype1, -i \
        infile2,unc2,etype2 or include argument: -d slab2database_location.'
        sys.exit(0)

    catalogs = []
    file_no = 1
    #Writes each file in filelist (events in bounds within each file) to output file
    for file in filelist:
        writetofile(file.filename, slab_, file.event_type, file.unc, args, catalogs,file_no)
        file_no = file_no+2

    # Makes rough plot of data for the region
    s2d_fnctns.slabplotter(args)

#Help/description and command line argument parser
if __name__=='__main__':
    desc = '''
        Writes a file in csv format with the information pertinent to Slab 2.0. Columns
            specified as:
        
        (lat,lon,depth,uncertainty,event-type,mag,time,P-azimuth,P-plunge,T-azimuth,
            T-plunge,strike1,dip1,rake1,strike2,dip2,rake2,slabStrike,slabDip,ID)
        
        Fields are represented by 'nan' where information is not available. Minimum input information
            includes columns labeled as lat,lon,depth. If the type of event is an earthquake, time 
            and mag are also required data columns.
        Where applicable, CMT information should be represented in tensorial form with 
            columns labeled: mrr,mtt,mpp,mrt,mrp and mtp.
        
        Expected event type input:
            Earthquake: EQ
            Receiver Function: RF
            Tomography: TO
            Active Source: AS
            Bathymetry: BA
            Earthquake Relocation: ER
            Centroid Moment Tensor: MT
            
        
        EXAMPLE: to compile slab and moment tensor information around Southern Peru in 2013 from
            original data files:
        
        s2d.py -b -85 -60 -25 -15 -s 2013-01-01 -e 2014-01-01 -f speru13slab.csv 
            -i isc-gem.csv,15,EQ -i csn_cat2.txt,15,EQ -i so_peru_rf.txt,10,RF 
            -i s_peru_to.tsv,10,TO
        
        EXAMPLE: to compile slab and moment tensor information around Southern Peru from all 
            available files in slab2database:
        
        s2d.py -b -85 -60 -25 -15 -d /some/directory/slab2database -f speru_slab.csv
        
        
        The database and original files can be added in the same call of s2d.py.
        
        Note that output file and lat/lon bounds are required arguments. If mag and time bounds 
            are not specified, the values are set to all magnitudes and dates from 1900 to present.
        
        Note that when specifying a search box that crosses the -180/180 meridian, specify 
            longitudes as you would if you were not crossing that meridian (i.e., lonmin=179, 
            lonmax=-179).  The program will resolve the discrepancy.
        
        If more than one earthquake dataset is provided, the files will be compared and matching 
            events will not be written to the slab file. The program prioritizes catalogs in 
            the order in which they are entered, writing the matching events from the first entry 
            and disregarding the associated match(es) in later entries.
        If more than one match is found for a single event, the program selects the closest match 
            and determines the remaining potential matches as independent.
        
        Where CMT information is provided, non-thrust earthquakes at depths shallower than 60 km 
            are filtered out.
        
        The file is saved, and more information can be appended to it by calling to s2d.py.
        
        A local copy of s2d_fnctns.py is required to run this program.
        
        '''
    
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    #required arguments
    parser.add_argument('-f','--output-file', dest='outFile', type=str, required=True,
                        help='Designated output file where data is written.')
    parser.add_argument('-b','--bounds', metavar=('lonmin','lonmax','latmin','latmax'),
                        dest='bounds', type=float, nargs=4, required=True,
                        help='Bounds to constrain event search [lonmin lonmax latmin latmax]')
    #optional arguments
    parser.add_argument('-i','--inputFiles', dest='input', type=infile, action='append',
                        help = 'List of input files with their associated uncertainty and event-type. \
                        -i input1,unc1,etype1 -i input2,unc2,etype2')
    parser.add_argument('-d','--database', dest='database', type=str,
                       help = 'Directory where all current slab 2.0 files are located')
    parser.add_argument('-s','--start-time', dest='startTime', type=maketime,
                        help='Start time for search (defaults to 01/01/1900).  YYYY-mm-dd or \
                        YYYY-mm-ddTHH:MM:SS')
    parser.add_argument('-e','--end-time', dest='endTime', type=maketime,
                        help='End time for search (defaults to current date/time).  YYYY-mm-dd or \
                        YYYY-mm-ddTHH:MM:SS')
    parser.add_argument('-m','--mag-range', metavar=('minmag','maxmag'),dest='magRange',
                        type=float,nargs=2,
                        help='Min/max (authoritative) magnitude to restrict search.')

                        
    pargs = parser.parse_args()
                        
                        
    main(pargs)
