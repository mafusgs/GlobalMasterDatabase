#pylint: disable=C0103

#necessary imports and modules
import csv
from datetime import datetime
#from time import strftime
from ast import literal_eval as make_tuple
from libcomcat.comcat import getEventData
#from collections import OrderedDict
#record start time for this script

t0 = datetime.utcnow()
t0_secs = (t0-datetime(1970, 1, 1)).total_seconds()

#select time interval of interest
#datetime(year,month,day)
start = datetime(1973, 1, 1)
finish = datetime.utcnow()

#define magnitude range of earthquakes to search over for shallow earthquakes
#(these are whatever magnitude range is defined in the catalogues)
min_mag = 5
max_mag = 9.9
magrange = (min_mag, max_mag)

#shallow earthquakes near surface may or may not be related to subduction
#request moment tensor components to find thrust faults, more likely SZ related
#depths in kilometers
min_sh = 0
max_sh = 70
depthrange_sh = (min_sh, max_sh)

#define deep depth range
min_dp = 70
max_dp = 750
depthrange_dp = (min_dp, max_dp)

#set grid spacing depending upon what is needed to run code quickly/smoothly
#an integer, in degrees, and must be small because SW Alaska has many quakes
grdspace = 3

#define a grid which separates data queries into reasonable chunks
#bounds should be (lonmin,lonmax,latmin,latmax)
#bounds are set for global search
#(-180,180) is max longitude, (-85,85) is max latitude defined by ComCat
ymax = 85
ymin = ymax - grdspace
xmin = -180
xmax = xmin + grdspace
range_degrees_lon = 180
total_degrees_lon = 2*range_degrees_lon
range_degrees_lat = 85
total_degrees_lat = 2*range_degrees_lat

#boxes describes the total number of grids you have chosen to divide the world into
boxes = (total_degrees_lat/grdspace)*(total_degrees_lon/grdspace)

#maxiter represents max number of iterations in the y direction (longitude direction)
maxiter = total_degrees_lon/grdspace

#count keeps track of iterations of longitude
#holds latmin/latmax steady while lonmin/lonmax changes across
#when max iterations in longitude have completed (gone across the globe)
#the latmin/latmix will adjust and lonmin/lonmax will also be reset.
#This process will continue until the number of boxes has been reached.
count = 0

#initialize a file to store the bounds to
#'w' implies writing to a new file
gridfile = open('boundaries.csv', 'w')
writer = csv.writer(gridfile)

for i in range(boxes):

    if count == maxiter-1:
        lonmax = xmax + grdspace*count
        lonmin = xmin + grdspace*count
        count = 0
        latmax = ymax
        latmin = ymin
        boundaries = (lonmin, lonmax, latmin, latmax)
        ymax = ymax - grdspace
        ymin = ymin - grdspace
        writer.writerow((boundaries))

    else:
        lonmax = xmax + grdspace*count
        lonmin = xmin + grdspace*count
        count = count+1
        latmax = ymax
        latmin = ymin
        boundaries = (lonmin, lonmax, latmin, latmax)
        writer.writerow((boundaries))

gridfile.close()

#keep track of the first time data is written to the files
line_no_sh = 0
line_no_dp = 0

with open('boundaries.csv') as boundsfile:

    for line in boundsfile:

        bounds = make_tuple(line)

        #to follow along with the progress of data querying
        print bounds

        #use getEventData from comcat.py to search for data of interest
        #define mags for shallow search to get moment tensor info
        #getComponents=True for MT (i.e. mrt) where available
        shallowlist, magmax1 = getEventData(bounds=bounds, starttime=start, endtime=finish, depthrange=depthrange_sh, magrange=magrange, getComponents=True)

        #add today's date to filename
        filename1part1 = 'cc_shallowquakes_'
        filename1part2 = datetime.now().strftime("%Y-%m-%d")
        filename1part3 = '_.csv'
        myfilename1 = filename1part1 + filename1part2 + filename1part3

        #'a' writes and appends to file if it already exists
        myfile1 = open(myfilename1, 'a')
        writer = csv.writer(myfile1)

        #labels rows with useful headers but only do this once
        if line_no_sh == 0:

            writer.writerow(('id_no', 'time', 'lat', 'lon', 'depth', 'mag', 'event_type', 'mrr', 'mtt', 'mpp', 'mrt', 'mrp', 'mtp', 'type', 'moment_lat', 'moment_lon', 'moment_depth', 'moment_duration'))
            line_no_sh = 1

            for i in range(len(shallowlist)):
                col = shallowlist[i]
                #define positions of each variable
                id_no = col['id'][0]
                time = col['time'][0]
                lat = col['lat'][0]
                lon = col['lon'][0]
                depth = col['depth'][0]
                mag = col['mag'][0]
                event_type = col['event-type'][0]
                mrr = col['mrr'][0]
                mtt = col['mtt'][0]
                mpp = col['mpp'][0]
                mrt = col['mrt'][0]
                mrp = col['mrp'][0]
                mtp = col['mtp'][0]
                type_mt = col['type'][0]
                moment_lat = col['moment-lat'][0]
                moment_lon = col['moment-lon'][0]
                moment_depth = col['moment-depth'][0]
                moment_duration = col['moment-duration'][0]

                print "there are shallow earthquakes here"

                #write this information to csv file
                writer.writerow((id_no, time, lat, lon, depth, mag, event_type, mrr, mtt, mpp, mrt, mrp, mtp, type_mt, moment_lat, moment_lon, moment_depth, moment_duration))

                #now perform a query for deep earthquakes
                #but do not limit by magnitude here
                #all types associated with the slab at depths exceeding 70 km

		#use getEventData again
                deeplist, magmax2 = getEventData(bounds=bounds, starttime=start, endtime=finish, depthrange=depthrange_dp)

                #add today's date to the file name
                filename2part1 = 'cc_deepquakes_'
                filename2part2 = datetime.now().strftime("%Y-%m-%d")
                filename2part3 = '_.csv'
                myfilename2 = filename2part1 + filename2part2 + filename2part3

                myfile2 = open(myfilename2, 'a')
                writer2 = csv.writer(myfile2)

                if line_no_dp == 0:

                    writer2.writerow(('id_no', 'time', 'lat', 'lon', 'depth', 'mag', 'event_type'))
                    line_no_dp = 1

                for i in range(len(deeplist)):
                    col = deeplist[i]
                    #define positions of each variable
                    id_no = col['id'][0]
                    time = col['time'][0]
                    lat = col['lat'][0]
                    lon = col['lon'][0]
                    depth = col['depth'][0]
                    mag = col['mag'][0]
                    event_type = col['event-type'][0]

                    print "there are deep earthquakes here"

                    #write this information to csv file
                    writer2.writerow((id_no, time, lat, lon, depth, mag, event_type))

#close the files
myfile1.close()
myfile2.close()

#sum the time to run this script
t1 = datetime.utcnow()
t1_secs = (t1-datetime(1970, 1, 1)).total_seconds()
total = t1_secs - t0_secs
print 'total time elapsed:', total, 'seconds'
