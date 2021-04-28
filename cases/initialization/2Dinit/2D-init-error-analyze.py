"""

TODO: Description.

"""

from geomtransport_error import *
import sys

def main():

    # Process mesh intersection output.  

    intersectType = sys.argv[1]

    templateCase = intersectType + "Domain"
    dataFileName = intersectType + ".dat"

    if (intersectType == "cci"):
        columnNames = "Nt Nb Ni k l Ev Te".split(" ")
    elif (intersectType == "smci"):
        columnNames = "Nt Nb Ni Ev Te".split(" ")

    # Filter out directories that contain dataFile
    dataDirs = filter(lambda x : os.path.isdir(x) and dataFileName in os.listdir(x),
                      os.listdir(os.curdir))

    resultDataFrame = pd.DataFrame(columns=columnNames)
    i = 0
    for dataDir in dataDirs:
        dataFrame = pd.read_csv(os.path.join(os.curdir,dataDir,dataFileName),
                               header=None,names=columnNames)

        # Reduce the dataFrame to a series.
        dataSeries = pd.Series([dataFrame[x].mean() for x in columnNames],
                               index=columnNames)
        # Insert the dataSeries into the results dataFrame
        resultDataFrame.loc[i] = dataSeries
        i = i + 1

    resultDataFrame = resultDataFrame.sort_values('Nb')

    print(resultDataFrame)

    # TODO: Save raw data to HDF5 with metadata.
    # TODO: Save processed data to HDF5 with metadata.  

    #print(latexData)
    #
    #tableFileName = table_file_name(templateCase) 
    #
    #latexFile = open(tableFileName + ".tex", "w")
    #latexFile.write(latexData)

    #csvFile = open(tableFileName + ".csv", 'w')
    #csvFile.write(data.to_csv()) 

if __name__=="__main__":
    main()
