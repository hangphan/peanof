#!/usr/bin/env python3
from __future__ import division
'''
    PaEdiatric ANthropometric measurement Outlier Flagging pipeline
    @author Dr Hang Phan
    Clinical Informatics Research Unit (CIRU)
    NIHR BRC Data Science Cross-cutting 
    University of Southampton
    hang.phan@soton.ac.uk
    31/2/2020
'''
import pandas as pd
import numpy as np
from optparse import  OptionParser
from ggplot import * 
import matplotlib.pyplot as plt
from statsmodels.api import OLS
from sklearn import linear_model
import sys, os, logging
import logging.handlers

formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger('Log')
ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.setLevel(logging.INFO)

mapType={"CHILD_HEIGHT":"HEIGHT", "CHILD_WEIGHT":"WEIGHT", "CHILD_BMI":"BMI"}
mapTypeR = {"WEIGHT":"CHILD_WEIGHT", "HEIGHT": "CHILD_HEIGHT", "BMI":"CHILD_BMI", 'W': 'CHILD_WEIGHT', 'H':'CHILD_HEIGHT', 'WEIG': 'CHILD_WEIGHT', 'HEIG':'CHILD_HEIGHT'}
WHOThreshold ={'CHILD_HEIGHT': [-6,6], 'CHILD_WEIGHT': [-6,5], 'CHILD_BMI': [-5,5]}
AgeAdultHeight=18
heightChangeThres=1
LRCutoffSD = {'CHILD_HEIGHT': 2.0, 'CHILD_WEIGHT':2.9}
_baseDir = os.path.dirname(os.path.realpath(sys.argv[0]))


class childMeasurement(object):
    def __init__(self):
        self.cid = 'ChildID'
        self.gender = None
        self.dfA = pd.DataFrame(columns = ['ID', 'AGE', 'MEASURE_DATE', 'MEASURE_TYPE', 'MEASURE_VAL', 'SDS' ])
        self.dob = None
        self.stdVal = {'CHILD_HEIGHT': None, 'CHILD_WEIGHT':None}  #standard deviation scores of errors during linear regression of SDS values
        self.coef = {'CHILD_HEIGHT': None, 'CHILD_WEIGHT':None}    # coefficient of regression line
        self.intercept = {'CHILD_HEIGHT': None, 'CHILD_WEIGHT':None} # intercept of regression line
    def set_cid(self, cid) :
        self.cid=cid
    def get_nMeasurement(self, mType):
        return len(self.dfA[self.dfA.MEASURE_TYPE==mType])
    
    def get_firstAge(self, mType):
        if self.get_nMeasurement(mType) >0:
            return self.dfA[self.dfA.MEASURE_TYPE==mType].sort_values('AGE').AGE.tolist()[0]
        return None
    
    def get_lastAge(self, mType):
        if self.get_nMeasurement(mType) >0:
            return self.dfA[self.dfA.MEASURE_TYPE==mType].sort_values('AGE').AGE.tolist()[-1]
        return None
    def get_followUpTime(self, mType):
        start, end = self.get_firstAge(mType), self.get_lastAge(mType)
        if start!= None:
            return round(end - start,2)
    def getSummaryStats(self,mType):
        mapX = {'CHILD_HEIGHT': 'HEIG', 'CHILD_WEIGHT': 'WEIG'}
        prefix = mapX[mType]
        # df = pd.DataFrame([[1,1.23,'Hello']], columns=list('ABC'))   
        out = pd.DataFrame([[self.cid, self.get_nMeasurement(mType), self.get_firstAge(mType), self.get_lastAge(mType), self.get_followUpTime(mType),
                       self.stdVal[mType], self.coef[mType],self.intercept[mType]]],
                       columns = ['CID',prefix + '_N', prefix+'_FIRST_AGE', prefix + '_LAST_AGE', prefix + 'FU_LEN', 
                               prefix + '_STDVAL', prefix + '_COEFF', prefix + '_INTERCEPT'])
        
        return out
    
    def getAllSummaryStats(self):
        s1 = self.getSummaryStats('CHILD_HEIGHT')
        s2 = self.getSummaryStats('CHILD_WEIGHT')
        d = s1.merge(s2, on = 'CID', how = 'left')
        d = d.round({"HEIG_STDVAL":2, "HEIG_COEFF":2, "HEIG_INTERCEPT":2, "WEIG_STDVAL":2,  "WEIG_COEFF":2, "WEIG_INTERCEPT":2}) 
       
        return d
      
    def setGender(self, gender):
        if gender not in ['M', 'F']:
            logger.error("Invalid gender value, has to be M for male or F for female")
            sys.exit()
        self.gender = gender
        return
        
    def setDob(self, thisDate): #date format is "dd/mm/yyyy"
        self.dob = pd.to_datetime(thisDate, dayfirst=True)
        return
        
    def convertMeasureType(self, mType):
        if mType in ['CHILD_HEIGHT', 'CHILD_WEIGHT']:
            return mType
        if mType in mapTypeR:
            return mapTypeR[mType]
        else:
            logger.error('measureType is not one of H,W, HEIGHT, WEIGHT, HEIG, WEIG, CHILD_HEIGHT, CHILD_WEIGHT. Quit here')
            sys.exit()
        return
    

    def addMeasurementWithDate(self,  measureDates, measureType, measureVals):
        mType = self.convertMeasureType(measureType)
        if len(measureDates) != len(measureVals):
            logger.error('Number of measure dates and measure values do not match. Quit here')
            sys.exit()
            
        dft = pd.DataFrame({
            'ID' :pd.Series([self.cid]*len(measureDates)),
            'AGE':pd.Series([]*len(measureDates)),
            'MEASURE_DATE':pd.Series([pd.to_datetime(d, dayfirst=True) for d in measureDates]),
            'MEASURE_TYPE':pd.Series([mType]*len(measureDates)),
            'MEASURE_VAL' :pd.Series(measureVals)
        })
        dft['AGE'] = dft.apply(lambda x: self.ageCalculation(x.AGE, x.MEASURE_DATE), axis = 1)
        cols = ['ID', 'AGE', 'MEASURE_DATE', 'MEASURE_TYPE', 'MEASURE_VAL']
        self.dfA = self.dfA.append(dft, sort=False).sort_values('AGE').reset_index()[cols]
        return
    
    def addMeasurementWithAge(self, ages, measureType, measureVals):
        mType = self.convertMeasureType(measureType)
        if len(ages) != len(measureVals):
            logger.error('Number of ages and measure values do not match. Quit here')
            sys.exit()

        dft = pd.DataFrame({
            'ID' : [self.cid]*len(ages),
            'AGE':ages,
            'MEASURE_DATE': []*len(ages),
            'MEASURE_TYPE':[mType]*len(ages),
            'MEASURE_VAL' : measureVals
                            
        })
        cols = ['ID', 'AGE', 'MEASURE_DATE', 'MEASURE_TYPE', 'MEASURE_VAL']
        self.dfA = self.dfA.append(dft, sort=False).sort_values('AGE').reset_index()
        return

    def ageCalculation(self, age, dom):
        if pd.isnull(age):
            if pd.isnull(self.dob):
                logger.debug('Age is not provided. No date of birth to calculate age at measurement. Ignoring this record')
                return None
            return round((pd.to_datetime(dom, dayfirst=True) - self.dob).days/365.25, 2)
        else:
            return age
        
    def addMeasurementFromFile(self, childID, fn):
        dft = None
        if fn.endswith('csv'):
            dft = pd.read_csv(fn)
        if fn.endswith('xlsx'):
            dft = pd.read_excel(fn)

        dft.columns = ['ID', 'AGE', 'MEASURE_DATE', 'MEASURE_TYPE', 'MEASURE_VAL']
        dft = dft[dft.ID==childID]
        self.addMeasurementFromDf(childID, dft)
        return 

    def addMeasurementFromDf(self,  dfx):
        dft = dfx[dfx.CID==self.cid]
        dft['MEASURE_TYPE'] = dft['MEASURE_TYPE'].apply(lambda x: self.convertMeasureType(x))
        dft['MEASURE_DATE'] = pd.to_datetime(dft['MEASURE_DATE'], dayfirst=True)
        dft['AGE'] = dft.apply(lambda x: self.ageCalculation(x.AGE, x.MEASURE_DATE), axis = 1)
        self.dfA = dft.copy().sort_values(['MEASURE_TYPE', 'AGE']).reset_index()

        return
        
    def calculateSDSVal(self, childSDSC):
        if pd.isnull(self.gender):
            logger.error('Gender is not provided. Quit here')
            sys.exit()
        self.dfA['SDS'] = self.dfA.apply(lambda x: childSDSC.calculateSDSVal(x.MEASURE_TYPE, x.AGE, self.gender, x.MEASURE_VAL), axis=1)
        self.dfA['FILTER_FLAG'] = self.dfA.apply(lambda x: self.checkWHOThreshold(x.MEASURE_TYPE, x.SDS), axis=1)

    def plotGrowthChart(self,measureType, childSDSC, fon):
        dft = self.dfA[self.dfA.MEASURE_TYPE==measureType]
        ages = dft.AGE
        measureVals = dft.MEASURE_VAL
        childSDSC.plotGrowthChart( measureType, ages, self.gender, measureVals, fon, self.cid)
        return
    
    def checkWHOThreshold(self, measureType, SDSVal):
        if SDSVal < WHOThreshold[measureType][0] or SDSVal > WHOThreshold[measureType][1]:
            return 'WHO'
        else:
            return 'PLAUSIBLE'
        
    
    def printMeasurements(self):
        print (self.dfA.sort_values(['MEASURE_TYPE', 'AGE']))
    
    def getOutDf(self):
        return self.dfA.sort_values(['MEASURE_TYPE', 'AGE'])

    def weightCheck(self):
        '''
        For those that passed other filters, check if the subsequent value has gain or drop 
        a sliding window of 1 year to filter extreme changes in weight measurements:
        1. 1 day : +- 25%
        2. 3 months : +- 40%
        3. 1 year: +-50%
        '''
        dft = self.dfA[(self.dfA.MEASURE_TYPE == 'CHILD_WEIGHT')&(self.dfA.FILTER_FLAG!='WHO')]
        daysInYear=365.25
        oneDay, threeMonths, oneYear = 1.0/daysInYear, 0.25, 1
        idxList = dft.index.tolist()
        weightPass = [idxList[0]]
        for idx in idxList[1:]:
            elapseT = dft.loc[idx, 'AGE']  - dft.loc[weightPass[-1], 'AGE']
            absChangePct = abs(dft.loc[idx, 'MEASURE_VAL'] -  dft.loc[weightPass[-1], 'MEASURE_VAL'])/dft.loc[weightPass[-1], 'MEASURE_VAL']*100
            if elapseT <= oneDay:
                if absChangePct >=25:
                    self.dfA.loc[idx, 'FILTER_FLAG'] = 'WEIG_ONEDAY'
                    pass
                else:
                    weightPass.append(idx)
                    continue
            elif elapseT <= threeMonths:
                if absChangePct >=40:
                    self.dfA.loc[idx, 'FILTER_FLAG'] = 'WEIG_THREEMONTHS'
                    pass
                else:
                    weightPass.append(idx)
                    continue
            elif elapseT <= oneYear:
                if absChangePct >=50:
                    self.dfA.loc[idx, 'FILTER_FLAG'] = 'WEIG_ONEYEAR'
                    pass
                else:
                    weightPass.append(idx)
                continue
            else:
                weightPass.append(idx)
        return 
    

    def heightCheck(self):
        '''
        For those that passed other filters, check if the subsequent value has gain or drop 
        a sliding window of 3 months year to filter extreme changes in height measurements:
        3 months: >15%
        others: any height decrease
        #the vectors have been sorted by age
        '''
        daysInYear=365.25
        threeMonths =  0.25
        dft = self.dfA[(self.dfA.MEASURE_TYPE == 'CHILD_HEIGHT')&(self.dfA.FILTER_FLAG!='WHO')]
        idxList = dft.index.tolist()
        heightPass = [idxList[0]]
        for idx in idxList[1:]:
            elapseT = dft.loc[idx, 'AGE']  - dft.loc[heightPass[-1], 'AGE']
            changePct = (dft.loc[idx, 'MEASURE_VAL'] -  dft.loc[heightPass[-1], 'MEASURE_VAL'])/dft.loc[heightPass[-1], 'MEASURE_VAL']*100
            if elapseT <= threeMonths and changePct >=15:
                self.dfA.loc[idx, 'FILTER_FLAG'] = 'HEIGHT_INC15PCT'
                continue
            if dft.loc[idx, 'MEASURE_VAL'] -  dft.loc[heightPass[-1], 'MEASURE_VAL']< -1:
                self.dfA.loc[idx, 'FILTER_FLAG'] = 'HEIGHT_DEC1CM'
                pass
            else:
                heightPass.append(idx) 
                continue
        return
    
    def filterByAdultHeight(self):
        dft =  self.dfA[(self.dfA.MEASURE_TYPE == 'CHILD_HEIGHT')&(~self.dfA.FILTER_FLAG.isin(['WHO', 'OLS_OUTLIER', 'OLS_FEW_REMAIN']))]
        dft=dft[dft.AGE>=18]
        heights = dft.MEASURE_VAL
        if len(heights) ==0:
            return 
        heightAtAdultHeightEst = np.median(heights)
        for idx, row in dft.iterrows():
            if abs(row.MEASURE_VAL - heightAtAdultHeightEst)>heightChangeThres:
                self.dfA.loc[idx, 'FILTER_FLAG'] = 'HEIGHT_ADULT_OUTLIER'
        return

    def heightDecreaseCheck(self):
        dft =  self.dfA[(self.dfA.MEASURE_TYPE == 'CHILD_HEIGHT')&(~self.dfA.FILTER_FLAG.isin(['WHO', 'OLS_OUTLIER', 'HEIGHT_ADULT_OUTLIER']))]
        if len(dft)==0:
            return
        heights = dft.MEASURE_VAL.tolist()
        ages = dft.AGE.tolist()
        heightDecreaseThres=heightChangeThres
        lastHeight = heights[0]
        for idx in range(1,len(heights)):
            if lastHeight - heights[idx]>heightDecreaseThres:
                self.dfA.loc[dft.index[idx], 'FILTER_FLAG'] = 'HEIGHT_DEC_OUTLIER'
            else:
                lastHeight = heights[idx]
        return 

    def LRFunc(self, measureType):
        '''
        Linear regression using OLS 
        cut-off leverage: 3k/n
        cut-off for influence: 1
        cut-off for DFFITS 2*sqrt(k/n)
        cut-off for DFBETAS 2/sqrt(n)  where k=1
        '''
        dft =  self.dfA[(self.dfA.MEASURE_TYPE == measureType)&(self.dfA.FILTER_FLAG!='WHO')].copy()
        reg = linear_model.LinearRegression()
        regression = OLS(dft.MEASURE_VAL,dft.AGE).fit()
        infl = regression.get_influence()
        test = regression.outlier_test()

        k=1
        N = len(dft)
        dft['OLS_BONFPVAL'] = test['bonf(p)']
        dft['OLS_STUDENTRES']= test['student_resid']
        dft['OLS_INFLUENCE'] = infl.summary_frame().cooks_d
        dft['OLS_DFFITS'] = infl.summary_frame().dffits
        dft['OLS_DFB_AGE']=infl.summary_frame().dfb_AGE
        dft['N'] = [N] * N 

        coL, coI, coDf1, coDf2 = 3.0*k/N, 1, 2*(k/N)**0.5, 2/(N**0.5)
        dft1 = dft[(abs(dft['OLS_INFLUENCE'])<=coI)&(abs(dft['OLS_DFFITS'])<=coDf1) &(abs(dft['OLS_DFB_AGE'])<=coDf2) ]
        
        if len(dft1) <=2:
            for idx,row in dft.iterrows():
                self.dfA.loc[idx, 'FILTER_FLAG'] = 'OLS_FEW_REMAIN'
            return
        
        reg.fit(dft1[['AGE']], dft1['SDS'])
        dft['pred1'] = reg.predict(dft[['AGE']])
        dft['diff1'] = dft['SDS'] - dft['pred1']
        stdVal = dft[dft.index.isin(dft1.index)].diff1.std()
        dft['STD_FOLD'] = dft.diff1/stdVal

        self.stdVal[measureType]= stdVal
        self.coef[measureType] = reg.coef_[0]
        self.intercept[measureType] = reg.intercept_
        
        for idx, row in dft.iterrows():
            if abs(row.STD_FOLD) <= LRCutoffSD[measureType]:
                self.dfA.loc[idx, 'FILTER_FLAG'] = 'PLAUSIBLE'
            else:
                self.dfA.loc[idx, 'FILTER_FLAG'] = 'OLS_OUTLIER'

        return
    
    
    def outlierFlaggingFull(self, measureType, childSDSC = None, fon = None):
        dft = self.dfA[(self.dfA.MEASURE_TYPE==measureType)  & (~self.dfA.FILTER_FLAG.isin(['WHO']))]
        if len(dft) ==0:
            return
        if len(dft) ==1:
            self.dfA.loc[dft.index[0], 'FILTER_FLAG'] = 'PLAUSIBLE'
            if childSDSC:
                self.plotGrowthChart(measureType, childSDSC, fon)
            return
        if len(dft) in [2,3]:
            if measureType == 'CHILD_HEIGHT':
                self.heightCheck()
            if measureType == 'CHILD_WEIGHT':
                self.weightCheck()
            if childSDSC:
                self.plotGrowthChart(measureType, childSDSC, fon)
            return

        self.LRFunc(measureType)
        if measureType == 'CHILD_HEIGHT':
            self.filterByAdultHeight()
            self.heightDecreaseCheck()

        if fon!= None:
            if childSDSC ==None:
                logger.info('Child SDS class not provided')
                sys.exit()
            if not fon.endswith('png'):
                logger.error('Filename not ending with png. Quit here')
                sys.exit()
            self.plotGrowthChart(measureType, childSDSC, fon)
            self.plotSDS(measureType, fon.replace('.png', '_sds.png'))
        return

    def plotSDS(self, measureType, fon):
        dft = self.dfA[self.dfA.MEASURE_TYPE == measureType]
        
        figure = plt.figure(figsize = (10,8))
        ax = figure.add_subplot(1,1,1)
        ax.scatter(dft.AGE, dft.SDS, color='blue', label='SDSVals',s=30)
        ax.set_xlabel("age")
        ax.set_ylabel("SDS values {0} of {1}".format(measureType, self.cid))
        ax.set_ylim(-7,7)
        ax.set_title("{0} SDS values".format(self.cid))
        ax.set_xlim(0,20)
   
        a, b = self.coef[measureType], self.intercept[measureType]
        x = np.array([dft.AGE.min(), dft.AGE.max()])
        f = lambda t: a*t + b
        ax.plot(x, f(x), c = 'green')
        f1 = lambda t: a*t + b + LRCutoffSD[measureType] * self.stdVal[measureType]
        f2 = lambda t: a*t + b - LRCutoffSD[measureType] * self.stdVal[measureType]
        ax.plot(x, f1(x), color = 'red', linestyle='dashed')
        ax.plot(x, f2(x), color = 'red', linestyle='dashed')
        plt.savefig(fon)
        return
    
class childSDS(object):
    def __init__(self, opts):
        self.refGrowthFile = opts.refg
        self.dfg =pd.read_csv(self.refGrowthFile, header=0)
        self.sds = opts.sds
        return 
    def getLowerClosest(self, myArr, myNumber):
        lower = min([ i for i in myArr if i <= myNumber], key=lambda x:abs(x-myNumber))
        return lower

    def calculateSDSValD(self, measureType, dob, measureDate, gender, val):
        age = (pd.to_datetime(measureDate, dayfirst=True) - pd.to_datetime(dob, dayfirst=True)).days/365.25
        return self.calculateSDSVal(measureType, age, gender, val)

    def checkInfos(self, measureType, gender):
        if gender not in ['M', 'F']:
            logger.error('Provided gender not in dictionary, should be M for male and F for female')
            sys.exit()
        mType = measureType
        if mType not in ['CHILD_HEIGHT', 'CHILD_WEIGHT']:
            
            if measureType not in mapTypeR:
                logger.error('Provided measureType not in dictionary, should be HEIGHT, WEIGHT, HEIG, WEIG, H or W')
                sys.exit()
            mType = mapTypeR[measureType]
        return mType
    
    def calculateSDSVal(self, measureType, age, gender, measureVal):
        mType = self.checkInfos(measureType, gender)
        dfSlice=self.dfg[((self.dfg.MEASURE_TYPE == mType) & (self.dfg.GENDER==gender))] #GRT: grow reference table
        ageVecInGRT = dfSlice['AGE'].tolist()
        #age not in range for SDS calculation
        if age > max(ageVecInGRT) or age < min(ageVecInGRT):
            logger.debug('Age {} not in range for SDS calculation.'.format(age))
            return None
        ageForSDS = self.getLowerClosest(ageVecInGRT, age)
        dfSlice.set_index(['MEASURE_TYPE', 'GENDER', 'AGE']) #data frame for grow reference table
        val = dfSlice.loc[ dfSlice['AGE']== ageForSDS]
        M,L,S = float(val.M.iloc[0]), float(val.L.iloc[0]), float(val.S.iloc[0])
        SDSVal = ((measureVal/M)**L-1)/(L*S)
        if self.sds:
            logger.debug("SDS value of the child {} measurement of {} is {}".format(measureType, measureVal, round(SDSVal,2)))
        return round(SDSVal,2)

    

    def plotGrowthChartD(self, measureType, dob, measureDates, gender, val, fon, cid = "childID"):
        dobD = pd.to_datetime(dob, dayfirst=True)
        ages = [(pd.to_datetime(d, dayfirst=True) - dobD).days/365.25 for d in measureDates]
        self.plotGrowthChart(measureType, ages, gender, val, fon)
        return
    
    def plotGrowthChart(self,  measureType, ages, gender, measureVals, fon, cid = "childID"): #plotThis could be output file name
        mType = self.checkInfos(measureType, gender)
        colsForVariable= ["CENTILE_0_4","CENTILE_2","CENTILE_9","CENTILE_25","CENTILE_50","CENTILE_75","CENTILE_91","CENTILE_98","CENTILE_99_6"]   
        dfS = self.dfg[((self.dfg.MEASURE_TYPE ==mType ) & (self.dfg.GENDER==gender))]
        dft = pd.DataFrame({'AGE': ages, mType: measureVals})
        dfS = dfS.append(dft, sort=False)
        if mType == 'CHILD_HEIGHT':
            ylabel = 'height (cm)'
        else:
            ylabel = 'weight (kg)'
        p = ggplot(aes(x='AGE'), data=dfS) + ylab(ylabel)
        p += ggtitle('Growth chart of child {}, {}, gender {}'.format(cid, mType, gender))
        for col in colsForVariable:
            p += geom_line(aes(y = col), data = dfS)
        p += geom_point(aes(y = mType, color = "red"))
        p.save(fon , width = 10, height=8, dpi = 300)
        return


class PEANOF(object):
    def __init__(self, sdsModule, fin, fon,  ids=None):
        self.fin = fin
        self.sdsModule = sdsModule
        self.fon=fon
        self.ids = None
        self.children = []
        self.dfA = None # anthropometric dataframe
        self.dfAO = None
        self.dfSO = None
        if ids:
            self.ids = ids

        if self.fon == None:
            self.fon= self.fin.replace('.csv', '_processed.csv').replace('.xlsx', '_processed.csv')
        if self.fon.endswith('csv'):
            pass
        else:
            self.fon += '.csv'
            
        if self.fin.endswith('xlsx'):
            self.dfA = pd.read_excel(self.fin)
        elif self.fin.endswith('csv'):
            self.dfA = pd.read_csv(self.fin)
        else:
            logger.error('Input file is not csv or excel (.xlsx) file. Please use the correct file format')
            sys.exit()
        self.readInputfiles()
        
    def readInputfiles(self):
        if len(self.dfA.columns)!=7:
            logger.error('Input file does not have the correct number of columns. The columns expected are: CID (ID of children), Date of birth, Gender, Age at measurement, Measure date, Measure type, Measure value. If the age at measurement is empty then the date of birth should be available.')
            sys.exit()
        self.dfA.columns = ['CID', 'DOB', 'GENDER', 'AGE', 'MEASURE_DATE', 'MEASURE_TYPE', 'MEASURE_VAL']
        self.demoDf = self.dfA[['CID', 'DOB', 'GENDER']].drop_duplicates().reset_index()
        self.dfA = self.dfA[['CID', 'AGE', 'MEASURE_DATE', 'MEASURE_TYPE', 'MEASURE_VAL']].drop_duplicates().reset_index()
            
            
    def makeChildrenList(self):
        if self.ids == None:
            self.ids = self.demoDf.CID.tolist()
        logger.info('Preparing data processing for {} individuals.'.format(len(self.ids)))
        for cid in self.ids:
            dft = self.dfA[self.dfA.CID == cid]
            dfdt = self.demoDf[self.demoDf.CID == cid]
            childModule = childMeasurement()
            childModule.set_cid(cid)
            childModule.setGender(dfdt.GENDER.tolist()[0])
            childModule.setDob(dfdt.DOB.tolist()[0])
            childModule.addMeasurementFromDf(dft)
            self.children.append(childModule)
        return
    
    def runOutlierDetection(self, fonPrefix=None):
        prefix = fonPrefix
        self.makeChildrenList()
        
        if len(self.ids)>5: #do not put plots to file if the number of children to process exceed 5
            prefix = None
        if prefix != None:
            logger.info('Processing and generating growth charts for {} individuals'.format(len(self.ids)))
        for idx, child in zip(self.ids, self.children):
            child.calculateSDSVal(self.sdsModule)
            if prefix:
                fon = prefix + '_{}_height.png'.format(idx)
                child.outlierFlaggingFull('CHILD_HEIGHT', self.sdsModule, fon)
                fon1 = prefix + '_{}_weight.png'.format(idx)
                child.outlierFlaggingFull('CHILD_WEIGHT', self.sdsModule, fon1)
            else:
                child.outlierFlaggingFull('CHILD_HEIGHT')
                child.outlierFlaggingFull('CHILD_WEIGHT')
        
        self.dfAO = pd.concat([child.getOutDf() for child in self.children])
        self.dfSO = pd.concat([child.getAllSummaryStats() for child in self.children])
        cols=['CID','AGE','MEASURE_DATE','MEASURE_TYPE','MEASURE_VAL','SDS','FILTER_FLAG']
        self.dfAO[cols].to_csv(self.fon, index=False)
        self.dfSO.to_csv(self.fon.replace('.csv', '_summary.csv'), index=False)
        logger.info('Output files: measurements with flags {}, and summary statistics per individual {}'.format(self.fon, self.fon.replace('.csv', '_summary.csv')))
        logger.info('Outlier flagging completed.')

if __name__ == "__main__":
    usage = "usage: python %prog [options] \n" 
    version = "%prog 0.1"
    
    parser = OptionParser(usage=usage, version=version)
    parser.add_option("-r", "--refg", dest="refg", type="string", default= _baseDir +'/'+ 'Growth_Reference_Data.csv', help="Growth reference table used for calculating SDS values from age, gender and measurement")
    parser.add_option("-n", "--number", dest="sds", type = 'int', default= 0, help = '\n0: calculate SDS values only for one individual, require age or dob, measuredates, and measurements  \n1: outlier flagging for the whole input file, require input file with ID, DOB, GENDER, AGE, MEASURE DATE, MEASURE TYPE and MEASURE VALUE. If no AGE information is available, both MEASURE DATE and DOB must be present. These are required to calculate SDS values. If --ids (-i) is set to a list of IDs then would only process the specified IDs in the dataset')
    parser.add_option("-f", "--fn", dest="fn", type="string", default=None, help="Name of file containing height and weight measurements with age and gender, can be xlsx file or csv file, and the columns must be in order of ID, DOB, GENDER, AGE, MEASURE DATE, MEASURE TYPE, MEASURE VALUE")
    parser.add_option("-o", "--fon", dest = "fon", type="string", default = None, help="Name of output file")
    parser.add_option("-i", "--ids", dest = "ids", type="string", default = None, help="List of IDs of individuals to process, comma separated list, if no setting, the program will process all")
    parser.add_option("-p", "--prefix", dest = "prefix", type="string", default = 'peanof', help="Prefix of output images for --ids option in processing a limited number of children")
    parser.add_option("-a", "--age", dest = "age", type="string", default = '', help="Age at measurement. Can have multiple age value separated by comma")
    parser.add_option("-m", "--measurement", dest = "measurement", type="string", default = '', help = "Measured value. Can have multiple measurement separated by comma")
    parser.add_option("-t", "--measureType", dest = "measureType", type = "string", default = "WEIGHT", help = "Type of measurement, WEIGHT in kg or HEIGHT in cm, default is WEIGHT")
    parser.add_option("-g", "--gender", dest = "gender", type = "string", default = None, help = "Gender (M for Male and F for Female)")
    parser.add_option("-d", "--dob", dest = "dob", type = "string", default = None, help = "Date of birth of the child (dd/mm/yyyy)")
    parser.add_option("-e", "--dom", dest = "dom", type = "string", default = '', help = "Date of measurement of the child (dd/mm/yyyy), can have multiple date of measurements separated by comma")

    (opts, args) = parser.parse_args()
    if opts.age:
        ages = list(map(float, opts.age.split(',')))
        vals = list(map(float, opts.measurement.split(',')))
        if len(ages) != len(vals):
            logger.error('Number of measurements and number of ages at measurement do not match. Quit here.')
            sys.exit()
    if opts.dom:
        doms = opts.dom.split(',')
        vals = list(map(float, opts.measurement.split(',')))
        if len(doms) != len(vals):
            logger.error('Number of measurements and number of dates do not match. Quit here.')
            sys.exit()
    
    sdsModule = childSDS(opts)        
    if opts.sds == 0:
        if opts.age:
            for (age, val) in zip(ages, vals):
                sdsModule.calculateSDSVal(opts.measureType, age, opts.gender, val)
            sdsModule.plotGrowthChart(opts.measureType,  ages, opts.gender, vals, opts.fon)
        else:
            for (dom, val) in zip(doms, vals):
                sdsModule.calculateSDSValD(opts.measureType, opts.dob, dom, opts.gender, val)
            sdsModule.plotGrowthChartD(opts.measureType, opts.dob, doms, opts.gender, vals, opts.fon)

    elif opts.sds ==1:
        if opts.ids:
            outlierFlagC = PEANOF(sdsModule, opts.fn, opts.fon, ids= opts.ids.split(','))
            outlierFlagC.runOutlierDetection(fonPrefix=opts.prefix)
        else:
            outlierFlagC = PEANOF(sdsModule, opts.fn, opts.fon)
            outlierFlagC.runOutlierDetection()                        


