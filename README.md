## peanof - PaEdiatric ANthropometric measurement Outlier Flagging pipeline

**Automated data cleaning of paediatric anthropometric data from longitudinal electronic health records: protocol and application to a large patient cohort**

*Hang T.T. Phan, Florina Borca, David Cable, James Batchelor, Justin H. Davies, Sarah Ennis*
*University of Southampton, Southampton, UK*

### Dependencies for direct use of the peanof.py script
The pipeline requires Python 3.7 and the following packages to run.

1. ggplot
2. sklearn
3. statsmodels
4. matplotlib

### Usages
    Usage: python peanof.py [options] 

    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      -r REFG, --refg=REFG  Growth reference table used for calculating SDS values
                        from age, gender and measurement
      -n SDS, --number=SDS   0: calculate SDS values only for one individual,
                            require age or dob, measuredates, and measurements
                            1: outlier flagging for the whole input file, require
                            input file with ID, DOB, GENDER, AGE, MEASURE DATE,
                            MEASURE TYPE and MEASURE VALUE. If no AGE information
                            is available, both MEASURE DATE and DOB must be
                            present. These are required to calculate SDS values.
                            If --ids (-i) is set to a list of IDs then would only
                            process the specified IDs in the dataset
      -f FN, --fn=FN        Name of file containing height and weight measurements
                            with age and gender, can be xlsx file or csv file, and
                            the columns must be in order of ID, DOB, GENDER, AGE,
                            MEASURE DATE, MEASURE TYPE, MEASURE VALUE
      -o FON, --fon=FON     Name of output file
      -i IDS, --ids=IDS     List of IDs of individuals to process, comma separated
                            list, if no setting, the program will process all
      -p PREFIX, --prefix=PREFIX
                            Prefix of output images for --ids option in processing
                            a limited number of children
      -a AGE, --age=AGE     Age at measurement. Can have multiple age value
                            separated by comma
      -m MEASUREMENT, --measurement=MEASUREMENT
                            Measured value. Can have multiple measurement
                            separated by comma
      -t MEASURETYPE, --measureType=MEASURETYPE
                            Type of measurement, WEIGHT in kg or HEIGHT in cm,
                            default is WEIGHT
      -g GENDER, --gender=GENDER
                            Gender (M for Male and F for Female)
      -d DOB, --dob=DOB     Date of birth of the child (dd/mm/yyyy)
      -e DOM, --dom=DOM     Date of measurement of the child (dd/mm/yyyy), can
                            have multiple date of measurements separated by comma
                        
### Example use cases
#### 1. Calculating SDS values for height or weight measurements

'python3 peanof.py -n 0 -d 1/01/2010 -e 12/3/2020,11/11/2019 -m 139,130 -g M -t HEIGHT   -o testNow.png'

#### 2. Outlier flagging of height/weight measurements for a few individuals from a big input file
The following will produce flagging results as well as the growth charts of 'Child1' and 'Child2' in the 'testData' folder

'python3 peanof.py -n 1 -f testData/test.csv -i Child1,Child2 -o testOut.csv -p testData/test'



#### 3. Outlier flagging of height/weight measurments for the whole dataset


### Docker usage
#### 1. Docker image build or pull
First download the Dockerfile to a folder of where you plan to build the docker image, then run the command:

'docker build -t peanof:1.0 .'

Or you can pull the peanof image from the public Docker repository 
'docker pull peanof:1.0'

#### 2. Running pipeline using docker image
##### a. Calculating SDS values input directly from commandline
    docker run --rm  \
      --name devtest \
      --mount source=`pwd`/testData,target=/data,type=bind \
      peanof:1.0 peanof.py -n 0 -d 1/01/2010 -e 12/3/2020,11/11/2019 -m 139,130 -g M -t HEIGHT   -o testNow.png

##### b. Running outlier flagging pipeline for an xlsx or csv inputfile containing height and weight measurements
    docker run --rm \
      --name devtest \
      --mount source=`pwd`/testData,target=/data,type=bind \
      peanof:1.0 peanof.py -n 1 -f test.csv