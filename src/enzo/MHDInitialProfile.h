/*
 * InitialProfile.h
 *
 *  Created on: Mar 18, 2019
 *      Author: cdev
 */

#ifndef SRC_ENZO_MHDINITIALPROFILE_H_
#define SRC_ENZO_MHDINITIALPROFILE_H_

#define LINE_MAX_LENGTH (512)
#define CMP_A_B(a, b) ((b>a)-(b<a))
#define SIGN_A(a) (CMP_A_B(0,a))

#define PROFILE_ASC_SORT (1)
#define PROFILE_DESC_SORT (-1)
#define PROFILE_UNSORTED (0)

#define PROFILE_EMPTY_LINE (0)
#define PROFILE_COMMENT_LINE (1)
#define PROFILE_COL_NAMES_LINE (2)
#define PROFILE_DATA_LINE (3)
#define PROFILE_TIME_LINE (4)

struct MHDInitialProfile
{
	long long nCols, nRows, nColsToKeep, nRowsAllocated, nRowsAllocateFirst, nRowsAllocateInc;
	char **colNames, **colNamesToKeep;
	long long *colSortingOrders, *colNumsToKeep;
	double requestedTime, frameTimeFound;
	double **colData;
	char* radiusColumnName;
	char* radialVelocityColumnName;
	char* densityColumnName;
	char* internalEnergyColumnName;
	char* temperatureColumnName;
	double* radiusData;
	double* radialVelocityData;
	double* densityData;
	double* temperatureData;
	double* internalEnergyData;
	double* internalEnergyRadiusData;
	long long nRowsInternalEnergy;
	double* pressureData;
	long long radiusIndex;
	long long radialVelocityIndex;
	long long densityIndex;
	long long temperatureIndex;
	long long internalEnergyIndex;
	long long radiusSortingOrder;

	/**
	 * Returns the position of the next whitespace character in s following start.
	 */
	char* strSkipSpace(char* s);
	/**
	 * Returns the position of the next non-whitespace character following start.
	 */
	char* strSkipNonSpace(char* s);
	/**
	 * Returns the count of non-empty words in a string separated by whitespace.
	 */
	long long countWords(char* s);
	/**
	 * Returns the index of a string in an array of strings with n elements
	 * or -1 if not found. NULLs in a are ignored.
	 */
	long long indexOfStr(char* s, char** a, long long n);
	/**
	 * Returns the index of the first occurrence of an integer
	 * in an array of integers or -1 if not found.
	 */
	long long indexOfInt(long long i, long long* a, long long n);
	/**
	 * Returns the would-be sorting order if nextElement was appended
	 * to the array prevElements with nPrevElements number of elements,
	 * which are already in prevSortingOrder.
	 */
	long long nextSortingOrder(double nextElement, double* prevElements, long long nPrevElements,
		long long prevSortingOrder);
	/**
	 * Returns the index of the first number in the xData array
	 * with nData elements, that is greater (smaller) or equal to x,
	 * in case xData is sorted in ascending (descending) order.
	 *
	 *  * Returns -1 if not found, i.e. x is greater (smaller) than
	 * all elements of xData.
	 *
	 * The sorting order is determined by xSortingOrder.
	 *
	 * For this application xData is assumed to be sorted without
	 * repetitions, i.e. strictly monotonic, even though this function
	 * may find a valid index in an non-monotonic or unsorted array
	 * for some search algorithms.
	 */
	long long findGTEIndex(double x, double xData[], long long nData, long long xSortingOrder);
	/**
	 * Interpolates y(x) using the data arrays xData and yData.
	 * If x is outside the range of xData, the value of
	 * the corresponding end element of yData is used
	 * as extrapolation.
	 *
	 * Uses findGTEIndex to find the right index of the bracket for x.
	 * Returns a non-zero error code if findGTEIndex fails,
	 * otherwise returns 0.
	 */
	long long interpolate(double* y, double yData[], double x, double xData[], long long nData,
		long long xSortingOrder);
	/**
	 * Returns the type of the line being processed.
	 */
	/**
	 * Parses a string of whitespace-separated words.
	 * colNames is used to store pointers to the words.
	 * Each word is allocated space dynamically and terminated by '\0'.
	 * Returns the count.
	 */
	long long parseColumnNames(char* line, char** colNames);

	/**
	 * The destructor.
	 */
	void free();

	void addColToKeep(char* colName, char** thisColName);

	/**
	 * Profile::init()
	 *
	 * Initializes field members to reflect an empty profile state with no data
	 * nor columns. Used by constructors and destructor.
	 */
	void init();
	void init(char* radiusColumnName, char* densityColumnName, char* InternalEnergyColumnName,
		char* temperatureColumnName, char* RadialVelocityColumnName);
	/**
	 * Allocates new colDolata, colNames and colSortingOrders
	 * initialized with NULLs.
	 * Doesn't allocate any data rows. See also: allocateMoreRows(...).
	 */
	long long allocateCols(long long nCols);

	/**
	 * Increases the allocation of the data rows arrays by a given number,
	 * preserving the data that is already in the arrays.
	 */
	long long allocateMoreRows(long long nMoreRows);

	/**
	 * Returns the index of a column name if existing or -1.
	 */
	long long findColIndex(char* colName);
	/**
	 * Returns a pointer to the column data if index and name exist.
	 */
	double* findCol(long long colIndex);
	/**
	 * Returns a pointer to the column data if index and name exist.
	 */
	double* findCol(char* colName);

	/**
	 * Interpolates y(x) bsaed on the column indices.
	 */
	long long interpolate(double* y, long long yColIndex, double x, long long xColIndex);
	/**
	 * Interpolates y(x) based on column names.
	 */
	long long interpolate(double* y, char* yname, double x, char* xname);
	long long interpolateDensity(double* y, double x);
	long long interpolateInternalEnergy(double* y, double x);
	long long interpolateTemperature(double* y, double x);
	long long interpolateRadialVelocity(double* y, double x);

	/**
	 * Turns column names in colNames into NULL,
	 * if not present in the colNamesToKepp or
	 * the colNumsToKeep lists.
	 */
	long long removeCols();

	/**
	 *
	 */
	void identifyNamedCols();

	/**
	 * Parses a columns header line to allocate the correct number of columns
	 * and save the column names. Returns the number of columns retained,
	 * which should be equal to nColsToKeep.
	 */
	long long processColNamesLine(char* line);
	/**
	 * Parses a data line as a whitespace-separated list of floats.
	 * Stores the data in the data rows arrays and increments nDataRows.
	 * Increases the allocation for data rows if necessary.
	 * Columns with NULL names are skipped.
	 */
	long long processDataLine(char* line, long long* nLinesWithExtraCols);
	/**
	 * Sets to NULL column names in varNames if they are not in the colNamesToKeep
	 * list. If colNamesToKeep is not NULL the colNumsToKeep is ignored. Otherwise
	 * sets to NULL col names
	 * @param colNamesToKeep
	 * @param colNumsToKeep
	 * @param nColsToKeep
	 * @return
	 */
	long long lineType(char* line);

	long long read(char* filename, char* format, double atTime);

	/**
	 * class Profile_PAH01 : public Profile
	 * This is a format allowing for profiles at different times.
	 * A profile snapshot in time starts with a TIME_LINE
	 * having three (3) floats, the first being the time.
	 * It is followed by a column header line and then multiple
	 * data rows which are the actual PROFILE_LINEs.
	 * There is a second column header for the same time and then
	 * more PROFILE_LINEs. The second group of profile lines
	 * is ignored.
	 */

	/**
	 * Overridden to distinguish between a TIME_LINE and a PROFILE_LINE,
	 * both of which look like DATA_LINEs to the base method.
	 */
	long long lineTypePAH01(char* line);
	/**
	 * A factory method for the PAH01 format. Creates a ProfilePAH01
	 * from the data in a given file, taking the profiles at the first
	 * time that exceeds or is equal to the requested *time*.
	 */
	long long readPAH01(char* filename, double atTime);
	/**
	 * Overridden to distinguish between a TIME_LINE and a PROFILE_LINE,
	 * both of which look like DATA_LINEs to the base method.
	 */
	long long lineTypePAH02(char* line);
	long long readPAH02(char* filename, double atTime);
};

#endif /* SRC_ENZO_MHDINITIALPROFILE_H_ */
