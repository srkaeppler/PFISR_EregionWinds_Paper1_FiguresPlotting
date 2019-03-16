"""
CalculateMean.py

Author: srkaeppler
Email: skaeppl [at] clemson [dot] edu
Date: 03/20/2018


Purpose of this function is to calculate the mean dictionaries and have a class that
handles that data processing.  I will reuse this class in plotting software or in python notebooks
I just want something consistent.

10/30/2018 - including the plasma flow and average plasma flows

03/15/2019
This is the version of the code that is submitted as part of
"Incoherent Scatter Radar Observations of Auroral E-region Thermospheric Winds: The Role of Plasma Forcing near Summer Solstice"
This code is in ever way a research code that is presented "as is."
You are strongly encouraged to contact me if you have any questions about what I did.

The purpose of this code is to be an overarching class that parses the monthy ISR
files and calculates the weight mean, mean, and median.  

"""

import tables
import numpy
import os
import datetime
import glob
import sys
import copy
sys.path.append('/Users/srkaeppler/research/data/NSF_PFISR_Eregion_NeutralWinds/ProcessingCode')
from tools.configreader.ConfigReader import *

class CalculateMeanWinds:
    def __init__(self, configFile):
        # parse the config file into what I need to do the processing
        self.ConfigReader = ConfigReader()
        self.configFile = configFile
        if configFile:
            if os.path.isfile(configFile):
                self.config = self.ConfigReader.read(configFile)
                print self.config
                # self._check_config()
            else:
                raise Exception('Error: Check config file parameters')

            # parse config file at this point

        else:
            raise Exception('Error: No Config File Specified - exiting')
            sys.exit()

        return

    def readafile(self,fname):
        '''
        Function written by nicolls to read
        '''
        h5file=tables.open_file(fname)
        output={}
        for group in h5file.walk_groups("/"):
            output[group._v_pathname]={}
            for array in h5file.list_nodes(group, classname = 'Array'):
                output[group._v_pathname][array.name]=array.read()
        h5file.close()

        return output

    def WeightedMeanNaN(self,x, dx):
        xdx = x/dx
        # bevington 4.17
        weightedMean = numpy.nansum(x/(dx**2), axis=0)/numpy.nansum(1./dx**2, axis=0)
        # bevington 4.19
        weightedVar = 1./numpy.nansum((1./dx**2), axis=0)
        weightedStd = numpy.sqrt(weightedVar)
        n = [numpy.where(numpy.isfinite(xdx[:,i]))[0].shape[0] for i in range(xdx.shape[1])]
        # print n
        return weightedMean, weightedStd, n

    def CalculateMeanWindDict(self, DecimalHoursTimeGrid,dt,DecimalHours,inDict):

        """
        Transfer data into variables - least amount of work
        """
        ZonalWind = inDict['ZonalWind']
        MeridWind = inDict['MeridWind']
        errZonalWind = inDict['errZonalWind']
        errMeridWind = inDict['errMeridWind']
        ZonalE = inDict['ZonalE']
        MeridE = inDict['MeridE']
        ZonalFlow = inDict['ZonalFlow']
        MeridFlow = inDict['MeridFlow']
        errZonalFlow = inDict['errZonalFlow']
        errMeridFlow = inDict['errMeridFlow']
        NeMean = inDict['NeMean']
        AE = inDict['AE']
        KP = inDict['KP']
        AP = inDict['AP']
        F107 = inDict['F107']
        MeanAltitude = inDict['MeanAltitude']

        # # added on 10/30/2018
        # ZonalFlowFregion = inDict['ZonalFlowFregion']
        # MeridFlowFregion = inDict['MeridFlowFregion']



        outDict = dict()


        # index 0 will be the weightedMean, indix 1 normal mean, index2 median
        labels = ['WeightedMean', 'NaNmean', 'NaNmedian']
        MeanZonalWinds = numpy.zeros([3,DecimalHoursTimeGrid.shape[0],MeanAltitude.shape[0]])*numpy.nan
        MeanMeridWinds = numpy.zeros([3,DecimalHoursTimeGrid.shape[0],MeanAltitude.shape[0]])*numpy.nan
        StdZonalWinds = numpy.zeros([2,DecimalHoursTimeGrid.shape[0],MeanAltitude.shape[0]])*numpy.nan
        StdMeridWinds = numpy.zeros([2,DecimalHoursTimeGrid.shape[0],MeanAltitude.shape[0]])*numpy.nan
        NArrZonalWinds = numpy.zeros([DecimalHoursTimeGrid.shape[0],MeanAltitude.shape[0]])*numpy.nan
        NArrMeridWinds = numpy.zeros([DecimalHoursTimeGrid.shape[0],MeanAltitude.shape[0]])*numpy.nan
        ZonalWindsRaw = numpy.zeros([DecimalHoursTimeGrid.shape[0],MeanAltitude.shape[0], 1000])*numpy.nan
        MeridWindsRaw = numpy.zeros([DecimalHoursTimeGrid.shape[0],MeanAltitude.shape[0], 1000])*numpy.nan

        # include the flows
        # use the electric field to determine this
        MeanZonalFlow = numpy.zeros([3,DecimalHoursTimeGrid.shape[0]])*numpy.nan
        MeanMeridFlow = numpy.zeros([3,DecimalHoursTimeGrid.shape[0]])*numpy.nan
        StdZonalFlow = numpy.zeros([2,DecimalHoursTimeGrid.shape[0]])*numpy.nan
        StdMeridFlow = numpy.zeros([2,DecimalHoursTimeGrid.shape[0]])*numpy.nan
        NArrZonalFlow = numpy.zeros([DecimalHoursTimeGrid.shape[0]])*numpy.nan
        NArrMeridFlow = numpy.zeros([DecimalHoursTimeGrid.shape[0]])*numpy.nan
        ZonalFlowRaw = numpy.zeros([DecimalHoursTimeGrid.shape[0],MeanAltitude.shape[0], 1000])*numpy.nan
        MeridFlowRaw = numpy.zeros([DecimalHoursTimeGrid.shape[0],MeanAltitude.shape[0], 1000])*numpy.nan

        NeMeanRaw = numpy.zeros([DecimalHoursTimeGrid.shape[0],MeanAltitude.shape[0], 1000])*numpy.nan
        AERaw = numpy.zeros([DecimalHoursTimeGrid.shape[0], 1000])*numpy.nan
        APRaw = numpy.zeros([DecimalHoursTimeGrid.shape[0], 1000])*numpy.nan
        F107Raw = numpy.zeros([DecimalHoursTimeGrid.shape[0], 1000])*numpy.nan
        KPRaw = numpy.zeros([DecimalHoursTimeGrid.shape[0], 1000])*numpy.nan

        # # added on 10/30/2018
        # MeanZonalFlowFregion = numpy.zeros([3,DecimalHoursTimeGrid.shape[0]])*numpy.nan
        # MeanMeridFlowFregion = numpy.zeros([3,DecimalHoursTimeGrid.shape[0]])*numpy.nan

        MeanZonalE = numpy.zeros([3,DecimalHoursTimeGrid.shape[0]])*numpy.nan
        MeanMeridE = numpy.zeros([3,DecimalHoursTimeGrid.shape[0]])*numpy.nan


        for itime in range(0,len(DecimalHoursTimeGrid),1):
            stime = DecimalHoursTimeGrid[itime]-dt/2
            etime = DecimalHoursTimeGrid[itime]+dt/2
            # print 'stime, etime, DecHoursTimeGrid', stime,etime,DecimalHoursTimeGrid[itime]
            if stime<0:
                stime+=24.0
                qtime = numpy.where( (DecimalHours>stime) | (DecimalHours<etime) )[0]
            elif etime>24.0:
                etime-=24.0
                qtime = numpy.where( (DecimalHours>stime) | (DecimalHours<etime) )[0]
            else:
                qtime = numpy.where( (DecimalHours > stime) & (DecimalHours < etime) )[0]

            # print 'stime,etime',stime,etime
            # print 'decimalHours',DecimalHours[qtime]
            # print 'qtime shape', qtime.shape
            # print 'ZonalWind', numpy.nanmean(ZonalWind[qtime,:], axis=0)

            #zonal winds

            iend = ZonalWind[qtime,:].shape[0]
            ZonalWindsRaw[itime,:,0:iend] = ZonalWind[qtime,:].T
            iend = MeridWind[qtime,:].shape[0]
            MeridWindsRaw[itime,:,0:iend] = MeridWind[qtime,:].T
            tmpMean,tmpStd,n = self.WeightedMeanNaN(ZonalWind[qtime,:], errZonalWind[qtime,:])
            MeanZonalWinds[0,itime,:] = tmpMean
            MeanZonalWinds[1,itime,:] = numpy.nanmean(ZonalWind[qtime,:], axis=0)
            MeanZonalWinds[2,itime,:] = numpy.nanmedian(ZonalWind[qtime,:],axis=0)
            StdZonalWinds[0,itime,:] = tmpStd
            Narr = numpy.sum(numpy.isfinite(ZonalWind[qtime,:]), axis=0)
            NArrZonalWinds[itime,:] = Narr
            StdZonalWinds[1,itime,:] = numpy.nanstd(ZonalWind[qtime,:], axis=0)/numpy.sqrt(Narr)

            NeMeanRaw[itime,:,0:iend] = NeMean[qtime,:].T
            # print 'itime,qtime', itime, AE[qtime], qtime
            AERaw[itime,0:iend] = AE[qtime]
            APRaw[itime,0:iend] = AP[qtime]
            KPRaw[itime,0:iend] = KP[qtime]
            F107Raw[itime,0:iend] = F107[qtime]


            # zonal flows
            tmpZonalFlow = -MeridE/0.5*1e4 #tesla
            # tmpMean,tmpStd,n = WeightedMeanNaN(tmpZonalFlow[qtime], tmpZonalFlow[qtime]*0.1)
            # MeanZonalFlow[0,itime,:] = tmpMean
            # print numpy.nanmean(tmpZonalFlow[qtime])
            MeanZonalFlow[1,itime] = numpy.nanmean(tmpZonalFlow[qtime])
            MeanZonalFlow[2,itime] = numpy.nanmedian(tmpZonalFlow[qtime])
            # StdZonalFlow[0,itime,:] = tmpStd
            Narr = numpy.sum(numpy.isfinite(tmpZonalFlow[qtime]))
            NArrZonalFlow[itime] = Narr
            StdZonalFlow[1,itime] = numpy.nanstd(tmpZonalFlow[qtime])/numpy.sqrt(Narr)

            # MeanZonalFlowFregion[1,itime] = numpy.nanmean(ZonalFlowFregion[qtime])
            # MeanZonalFlowFregion[2,itime] = numpy.nanmedian(ZonalFlowFregion[qtime])


            # meridional winds
            tmpMean,tmpStd,n = self.WeightedMeanNaN(MeridWind[qtime,:], errMeridWind[qtime,:])
            MeanMeridWinds[0,itime,:] = tmpMean
            MeanMeridWinds[1,itime,:] = numpy.nanmean(MeridWind[qtime,:], axis=0)
            MeanMeridWinds[2,itime,:] = numpy.nanmedian(MeridWind[qtime,:],axis=0)
            StdMeridWinds[0,itime,:] = tmpStd
            Narr = numpy.sum(numpy.isfinite(MeridWind[qtime,:]), axis=0)
            NArrMeridWinds[itime,:] = Narr
            StdMeridWinds[1,itime,:] = numpy.nanstd(MeridWind[qtime,:], axis=0)/numpy.sqrt(Narr)

            # meridional flows
            tmpMeridFlow = ZonalE/0.5*1e4 # convert to Tesla
            # tmpMean,tmpStd,n = WeightedMeanNaN(tmpMeridFlow[qtime], tmpMeridFlow[qtime]*0.1)
            # MeanMeridFlow[0,itime,:] = tmpMean
            MeanMeridFlow[1,itime] = numpy.nanmean(tmpMeridFlow[qtime])
            MeanMeridFlow[2,itime] = numpy.nanmedian(tmpMeridFlow[qtime])
            # StdMeridFlow[0,itime,:] = tmpStd
            Narr = numpy.sum(numpy.isfinite(tmpMeridFlow[qtime]))
            NArrZonalFlow[itime] = Narr
            StdMeridFlow[1,itime] = numpy.nanstd(tmpMeridFlow[qtime])/numpy.sqrt(Narr)

            # MeanMeridFlowFregion[1,itime] = numpy.nanmean(MeridFlowFregion[qtime])
            # MeanMeridFlowFregion[2,itime] = numpy.nanmedian(MeridFlowFregion[qtime])

            # electric fields
            MeanZonalE[1,itime] = numpy.nanmean(ZonalE[qtime])
            MeanZonalE[2,itime] = numpy.nanmedian(ZonalE[qtime])
            MeanMeridE[1,itime] = numpy.nanmean(MeridE[qtime])
            MeanMeridE[2,itime] = numpy.nanmedian(MeridE[qtime])


        outDict['MeanMeridWinds'] = MeanMeridWinds
        outDict['StdMeridWinds'] = StdMeridWinds
        outDict['MeanZonalWinds'] = MeanZonalWinds
        outDict['StdZonalWinds'] = StdZonalWinds

        outDict['MeanMeridFlow'] = MeanMeridFlow
        outDict['StdMeridFlow'] = StdMeridFlow

        outDict['MeanZonalFlow'] = MeanZonalFlow
        outDict['StdZonalFlow'] = StdMeridFlow

        outDict['NarrMeridWinds'] = NArrMeridWinds
        outDict['NArrZonalWinds'] = NArrZonalWinds

        outDict['NarrMeridFlow'] = NArrMeridFlow
        outDict['NArrZonalFlow'] = NArrZonalFlow

        outDict['DecimalHoursTimeGrid'] = DecimalHoursTimeGrid
        outDict['MeanAltitude'] = MeanAltitude

        outDict['ZonalE'] = ZonalE
        outDict['MeridE'] = MeridE

        outDict['MeanZonalE'] = MeanZonalE
        outDict['MeanMeridE'] = MeanMeridE

        outDict['ZonalWindsRaw'] = ZonalWindsRaw
        outDict['MeridWindsRaw'] = MeridWindsRaw

        outDict['NeMeanRaw'] = NeMeanRaw
        outDict['APRaw'] = APRaw
        outDict['AERaw'] = AERaw
        outDict['KPRaw'] = KPRaw
        outDict['F107Raw'] = F107Raw

        # added on 10/30/2018
        # outDict['MeanZonalFlowFregion'] = MeanZonalFlowFregion
        # outDict['MeanMeridFlowFregion'] = MeanMeridFlowFregion

        # print '\n \n \n'
        # print 'NeMean'
        # print NeMeanRaw
        # print 'end NeMean'

        return outDict


    def main(self,fnames):

        """
        Parse a config file to get out the parameters we want
        """


        """
        Extract the data and calculate the means in a few different ways
        """


        ReturnDict = dict()
        with open('DataAttrition.txt', 'w') as f:
            # f.write("nTotal,nNan, nUseableInital, nAttrition, nUsableFinal \n ")
            for i in range(len(fnames)):
                InternalDict = dict()
                fname=fnames[i]
                print 'fname', fname
                output=self.readafile(os.path.join(fname))

                onameIndividualTime = fname.split('/')[-1][:-3]+'_IndividualTime.pdf'
                odir = '/'.join(fname.split('/')[0:-1])
                TimeStr = fname.split('/')[-1].split('_')[0]
                # outFile = os.path.join(odir,oname)
                # outFileContour = os.path.join(odir,onameContour)
                # outFilePolar = os.path.join(odir,onamePolar)

                # outfilefull = os.path.join(dpath,oname)
                # outfilefull2 = os.path.join(dpath, 'Median','%04d.png'%count)

                # print i
                # if i==0:
                Altitude=output['/Winds']['Altitude'];
                if Altitude.ndim == 2:
                    InternalDict['MeanAltitude']=numpy.mean(Altitude,axis=1)
                elif Altitude.ndim == 1:
                    InternalDict['MeanAltitude'] = Altitude



                InternalDict['ZonalWind'] = output['/Winds']['WindGeo'][:,:,0]
                InternalDict['MeridWind'] = output['/Winds']['WindGeo'][:,:,1]
                InternalDict['errZonalWind'] =output['/Winds']['errWindGeo'][:,:,0]
                InternalDict['errMeridWind'] =output['/Winds']['errWindGeo'][:,:,1]
                InternalDict['ZonalFlow'] = output['/VectorVels']['VestGmag'][:,:,0]
                InternalDict['MeridFlow'] = output['/VectorVels']['VestGmag'][:,:,1]
                InternalDict['errZonalFlow'] = output['/VectorVels']['errVestGmag'][:,:,0]
                InternalDict['errMeridFlow'] = output['/VectorVels']['errVestGmag'][:,:,1]
                InternalDict['time'] = output['/Time']['UnixTime'];
                InternalDict['mlt'] = output['/Time']['MLTDecHrs']
                InternalDict['slt'] = output['/Time']['LocalDecHrs']
                InternalDict['ZonalE'] = output['/ElectricFields']['Efield'][:,0]#.read()']
                InternalDict['MeridE'] = output['/ElectricFields']['Efield'][:,1]
                InternalDict['NeMean'] = output['/Ne']['MeanNeFitted']

                InternalDict['AE'] = output['/GeophysicalParameters']['AE']
                InternalDict['KP'] = output['/GeophysicalParameters']['KP']
                InternalDict['AP'] = output['/GeophysicalParameters']['AP']
                InternalDict['F107'] = output['/GeophysicalParameters']['F107']

                # InternalDict['ZonalFlowFregion'] = output['/Fregion']['VestGmag_300km'][:,0]
                # InternalDict['MeridFlowFregion'] = output['/Fregion']['VestGmag_300km'][:,1]

                # else:
                #     InternalDict['ZonalWind'] = numpy.concatenate((InternalDict['ZonalWind'],output['/Winds']['WindGeo'][:,:,0]),axis=0)
                #     InternalDict['MeridWind'] = numpy.concatenate((InternalDict['MeridWind'],output['/Winds']['WindGeo'][:,:,1]),axis=0)
                #     InternalDict['ZonalFlow'] = numpy.concatenate((InternalDict['ZonalFlow'],output['/VectorVels']['VestGmag'][:,:,0]),axis=0)
                #     InternalDict['MeridFlow'] = numpy.concatenate((InternalDict['MeridFlow'],output['/VectorVels']['VestGmag'][:,:,1]),axis=0)
                #     InternalDict['errZonalWind'] = numpy.concatenate((InternalDict['errZonalWind'],output['/Winds']['errWindGeo'][:,:,0]),axis=0)
                #     InternalDict['errMeridWind'] = numpy.concatenate((InternalDict['errMeridWind'],output['/Winds']['errWindGeo'][:,:,1]),axis=0)
                #     InternalDict['errZonalFlow'] = numpy.concatenate((InternalDict['errZonalFlow'], output['/VectorVels']['errVestGmag'][:,:,0]), axis=0)
                #     InternalDict['errMeridFlow'] = numpy.concatenate((InternalDict['errMeridFlow'], output['/VectorVels']['errVestGmag'][:,:,1]), axis=0)
                #     InternalDict['time'] = numpy.concatenate((InternalDict['time'],output['/Time']['UnixTime']),axis=0)
                #     InternalDict['mlt'] = numpy.concatenate((InternalDict['mlt'], output['/Time']['MLTDecHrs']),axis=0)
                #     InternalDict['slt'] = numpy.concatenate((InternalDict['slt'], output['/Time']['LocalDecHrs']),axis=0)
                #     InternalDict['ZonalE'] = numpy.concatenate((InternalDict['ZonalE'],output['/ElectricFields']['Efield'][:,0]), axis=0)#.read()']
                #     InternalDict['MeridE'] = numpy.concatenate((InternalDict['MeridE'],output['/ElectricFields']['Efield'][:,1]), axis=0)
                #     InternalDict['NeMean'] = numpy.concatenate((InternalDict['NeMean'],output['/Ne']['MeanNeFitted']), axis=0)
                #     InternalDict['AE'] = numpy.concatenate((output['/GeophysicalParameters']['AE'],output['/GeophysicalParameters']['AE']),axis=0)
                #     InternalDict['KP'] = numpy.concatenate((output['/GeophysicalParameters']['KP'], output['/GeophysicalParameters']['KP']), axis=0)
                #     InternalDict['AP'] = numpy.concatenate((output['/GeophysicalParameters']['AP'], output['/GeophysicalParameters']['AE']), axis=0)
                #     InternalDict['F107'] = numpy.concatenate((output['/GeophysicalParameters']['F107'], output['/GeophysicalParameters']['F107']), axis=0)

            # filter by velocity
            # Ibad=numpy.where( (InternalDict['ZonalWind'] > 500.0) | (InternalDict['ZonalWind'] < -500. ) | \
            #                   (InternalDict['MeridWind'] > 500.0) | (InternalDict['MeridWind'] < -500.) | \
            #                   (InternalDict['errZonalWind'] > 500.) | (InternalDict['errMeridWind'] > 500.) | \
            #                   (InternalDict['NeMean'] < 1.e11) \
            #                  )
                nTot = float(numpy.ravel(InternalDict['ZonalWind']).shape[0])
                print numpy.where(numpy.isfinite(numpy.ravel(InternalDict['ZonalWind'])) == False)[0]
                print 'nnan raw', numpy.where(numpy.isnan(numpy.ravel(InternalDict['ZonalWind']))==True)

                # note if you look at the raw output where outputs a tuple of two arrays
                # those two arrays combined will double count.

                nNan = float(numpy.where(numpy.isnan(numpy.ravel(InternalDict['ZonalWind'])) == True)[0].shape[0])
                qnan = numpy.where(numpy.isnan(numpy.ravel(InternalDict['ZonalWind'])) == True)
                print numpy.ravel(InternalDict['ZonalWind'])[qnan]

                qUse = numpy.where(numpy.isfinite(numpy.ravel(InternalDict['ZonalWind']))==True)
                print numpy.ravel(InternalDict['ZonalWind'])[qUse]

                nUsableInital = float(numpy.where(numpy.isfinite(numpy.ravel(InternalDict['ZonalWind']))==True)[0].shape[0])
                print 'Total starting Data,', numpy.ravel(InternalDict['ZonalWind']).shape[0], numpy.ravel(InternalDict['NeMean']).shape[0]
                # print 'Data initially Nans', numpy.ravel(numpy.where(numpy.isfinite(InternalDict['ZonalWind']) == False)).shape[0], numpy.ravel(numpy.where(numpy.isnan(InternalDict['ZonalWind']) == True)).shape[0]/nTot
                # print 'Usable Data', numpy.ravel(InternalDict['ZonalWind']).shape[0]-numpy.ravel(numpy.where(numpy.isnan(InternalDict['ZonalWind']) == True)).shape[0]
                print 'nUsableInital, nNan', nUsableInital, nNan
                nTotalCheck = 0.
                nUsableCheck = 0.
                f.write('%s \n'%TimeStr)
                f.write('Initial Usable Data \n')
                for i in range(InternalDict['ZonalWind'].shape[1]):
                    nUsableInitalAlt = float(numpy.where(numpy.isfinite(numpy.ravel(InternalDict['ZonalWind'][:,i]))==True)[0].shape[0])
                    print 'alt, initial usable data', InternalDict['MeanAltitude'][i],nUsableInitalAlt, InternalDict['ZonalWind'][:,i].shape
                    f.write('alt, %0.1f, %0.1f \n'%(InternalDict['MeanAltitude'][i],nUsableInitalAlt))
                    nTotalCheck = InternalDict['ZonalWind'][:,i].shape[0]+nTotalCheck
                    nUsableCheck = nUsableInitalAlt + nUsableCheck
                print 'initial nTotal check nUsable+nNan, nTot', nUsableInital+nNan, nTot
                print 'CHECK nUsableInital, nTotal', nUsableInital, nUsableCheck, nTot, nTotalCheck
                Ibad=numpy.where( (InternalDict['ZonalWind'] > self.config['DATADISCARD']['ZonalWindMax']) | \
                                  (InternalDict['ZonalWind'] < self.config['DATADISCARD']['ZonalWindMin'] ) | \
                                  (InternalDict['MeridWind'] > self.config['DATADISCARD']['ZonalWindMax']) | \
                                  (InternalDict['MeridWind'] < self.config['DATADISCARD']['ZonalWindMin']) | \
                                  (InternalDict['errZonalWind'] > self.config['DATADISCARD']['ErrorZonalWind'])  | \
                                  (InternalDict['errMeridWind'] > self.config['DATADISCARD']['ErrorMeridWind']) | \
                                  (InternalDict['NeMean'] < self.config['DATADISCARD']['NeMean']) \
                                 )



                # f.write("nTotal,nNan, nUseableInital, nAttrition, nUsableFinal \n ")

                InternalDict['ZonalWind'][Ibad]=numpy.nan
                InternalDict['MeridWind'][Ibad]=numpy.nan
                InternalDict['errZonalWind'][Ibad]=numpy.nan
                InternalDict['errMeridWind'][Ibad]=numpy.nan
                InternalDict['ZonalFlow'][Ibad] = numpy.nan
                InternalDict['MeridFlow'][Ibad] = numpy.nan
                InternalDict['errZonalFlow'][Ibad] = numpy.nan
                InternalDict['errMeridFlow'][Ibad] = numpy.nan
                InternalDict['NeMean'][Ibad] = numpy.nan

                # print 'Data After Filtering Nans', numpy.ravel(numpy.where(numpy.isnan(InternalDict['ZonalWind']) == True)).shape[0]
                # print 'Number of finite data points', numpy.ravel(numpy.where(numpy.isfinite(InternalDict['ZonalWind']) == True)).shape[0]
                # nUsableFinal = float(numpy.ravel(numpy.where(numpy.isfinite(InternalDict['ZonalWind']) == True)).shape[0])
                nUsableFinal = float(numpy.where(numpy.isfinite(numpy.ravel(InternalDict['ZonalWind']))==True)[0].shape[0])
                nNanFinal = float(numpy.where(numpy.isnan(numpy.ravel(InternalDict['ZonalWind'])) == True)[0].shape[0])
                print 'Ibad size', numpy.ravel(Ibad[0]).shape[0]
                print 'nUsableFinal, nNanFinal,nTot check', nUsableFinal, nNanFinal, nUsableFinal+nNanFinal, numpy.ravel(InternalDict['ZonalWind']).shape[0]

                print 'total ending data number', numpy.ravel(InternalDict['ZonalWind']).shape[0]
                nTotalCheck = 0.
                nUsableCheck = 0.
                f.write('\n Final Usable data \n')
                for i in range(InternalDict['ZonalWind'].shape[1]):
                    nUsableInitalAlt = float(numpy.where(numpy.isfinite(numpy.ravel(InternalDict['ZonalWind'][:,i]))==True)[0].shape[0])
                    print 'alt, initial usable data', InternalDict['MeanAltitude'][i],nUsableInitalAlt, InternalDict['ZonalWind'][:,i].shape
                    TotalCheck = InternalDict['ZonalWind'][:,i].shape[0]+nTotalCheck
                    nUsableCheck = nUsableInitalAlt + nUsableCheck
                    f.write('alt, %0.1f, %0.1f \n'%(InternalDict['MeanAltitude'][i],nUsableInitalAlt))
                print 'initial nTotal check nUsable+nNan, nTot', nUsableInital+nNan, nTot
                print 'CHECK nUsableInital, nUsableCheck, nTot, nTotalCheck ', nUsableInital, nUsableCheck, nTot, nTotalCheck
                print '\n \n'

                # do not need to do Ibad for AE,KP,F107, etc.
                f.write( '-------------Hours-------------------- \n')
                f.write('%s \n'%TimeStr)
                # try to determine how much data in terms of time
                for i in range(InternalDict['ZonalWind'].shape[1]):
                    qGoodData = numpy.where(numpy.isfinite(numpy.ravel(InternalDict['ZonalWind'][:,i]))==True)[0]
                    # print 'test WindData', numpy.ravel(InternalDict['ZonalWind'][:,i])[qGoodData]
                    print 'qgoodDataShape', qGoodData.shape
                    tmpHours = numpy.sum(InternalDict['time'][qGoodData,1]-InternalDict['time'][qGoodData,0])/3600.
                    totalTime = numpy.nansum(InternalDict['time'][:,1]-InternalDict['time'][:,0])/3600.
                    # print 'alt, initial usable data', InternalDict['MeanAltitude'][i],nUsableInitalAlt, InternalDict['ZonalWind'][:,i].shape
                    f.write('alt, %0.1f, %0.1f, %0.1f \n'%(InternalDict['MeanAltitude'][i], tmpHours,totalTime))
                f.write('\n ########################################## \n ')
                # in case I need to check that the filtering is working
                # for i in range(ZonalWind.shape[0]):
                #     print i, numpy.nanmax(ZonalWind[i,:]), numpy.nanmin(ZonalWind[i,:])


                MeanTime=numpy.nanmean(InternalDict['time'],axis=1)
                # filter out all times with nans
                qnan = numpy.where(numpy.isnan(MeanTime) == False)[0]

                # print 'MeanTime.shape', MeanTime.shape, mlt.shape, slt.shape

                if (MeanTime.shape[0] == InternalDict['mlt'].shape[0]) & (MeanTime.shape[0] == InternalDict['slt'].shape[0]):
                    MeanTime = MeanTime[qnan]
                    InternalDict['mlt'] = InternalDict['mlt'][qnan]
                    InternalDict['slt'] = InternalDict['slt'][qnan]
                else:
                    raise ValueError ("Wrong dimensions on time arrays")


                InternalDict['ZonalWind'] = InternalDict['ZonalWind'][qnan,:]
                InternalDict['MeridWind'] = InternalDict['MeridWind'][qnan,:]
                InternalDict['errZonalWind'] = InternalDict['errZonalWind'][qnan,:]
                InternalDict['errMeridWind'] = InternalDict['errMeridWind'][qnan,:]
                InternalDict['ZonalE'] = InternalDict['ZonalE'][qnan]
                InternalDict['MeridE'] = InternalDict['MeridE'][qnan]
                InternalDict['ZonalFlow'] = InternalDict['ZonalFlow'][qnan,:]
                InternalDict['MeridFlow'] = InternalDict['MeridFlow'][qnan,:]
                InternalDict['errZonalFlow'] = InternalDict['errZonalFlow'][qnan,:]
                InternalDict['errMeridFlow'] = InternalDict['errMeridFlow'][qnan,:]
                InternalDict['NeMean'] = InternalDict['NeMean'][qnan,:]

                # # added on 10/30/2018
                # InternalDict['ZonalFlowFregion'] = InternalDict['ZonalFlowFregion'][qnan]
                # InternalDict['MeridFlowFregion'] = InternalDict['MeridFlowFregion'][qnan]

                print 'AE', len(qnan), InternalDict['AE'].shape
                InternalDict['AE'] = InternalDict['AE'][qnan]
                InternalDict['KP'] = InternalDict['KP'][qnan]
                InternalDict['AP'] = InternalDict['AP'][qnan]
                InternalDict['F107'] = InternalDict['F107'][qnan]

                # have some sort of filtering
                # need to fiilter out wind estimates > 500 or 100 m/s

                # new function which will basically calculate the mean and then plot the data

                """
                Setting up the time grid
                """
                # dminute = self.config['TIME']['dMinutes']
                dhours = self.config['TIME']['TimeIntervalMinutes']/60.
                dt = self.config['TIME']['TimeIntervalLengthMinutes']/60.
                DecimalHoursTimeGrid =numpy.arange(0,24,dhours)
                DecimalTime = numpy.array([datetime.datetime.utcfromtimestamp(t) for t in MeanTime])
                DecimalHours = numpy.array([t.hour+t.minute/60.+t.second/3600. for t in DecimalTime])
                InternalDict['ut'] = DecimalHours




                # utTimeDict = CalculateMeanWindDict(DecimalHoursTimeGrid,dt,DecimalHours,MeanTime,\
                #                           ZonalWind,MeridWind,errZonalWind, \
                #                           errMeridWind,ZonalE, MeridE,MeanAltitude, \
                #                           ZonalFlow,MeridFlow,errZonalFlow,errMeridFlow,NeMean)

                outDict = dict()
                outDict['ut'] = self.CalculateMeanWindDict(DecimalHoursTimeGrid,dt,InternalDict['ut'],InternalDict)
                outDict['slt'] = self.CalculateMeanWindDict(DecimalHoursTimeGrid,dt,InternalDict['slt'],InternalDict)
                outDict['mlt'] = self.CalculateMeanWindDict(DecimalHoursTimeGrid,dt,InternalDict['mlt'],InternalDict)

                ReturnDict[TimeStr] = outDict
            f.close()
        # mltTimeDict = self.CalculateMeanWindDict(DecimalHoursTimeGrid,dt,mlt,MeanTime,\
        #                           ZonalWind,MeridWind,errZonalWind, \
        #                           errMeridWind,ZonalE, MeridE,MeanAltitude,\
        #                           ZonalFlow,MeridFlow,errZonalFlow,errMeridFlow,\
        #                           NeMean, AE, AP,KP, F107)
        return ReturnDict

if __name__=='__main__':
    print 'foo'
    configfile = '/Users/srkaeppler/research/data/NSF_PFISR_Eregion_NeutralWinds/MeanProfiles/meanprofile.ini'
    calc = CalculateMeanWinds(configfile)

    fnames = ['/Users/srkaeppler/Dropbox/research/data/NSF_PFISR_Eregion_NeutralWinds_SharedData/v0.4.2.2018.10.15/MonthlyWinds/201106_Winds.v0.4.2.2018.10.15.h5']
    testDict = calc.main(fnames)
    print testDict.keys()
