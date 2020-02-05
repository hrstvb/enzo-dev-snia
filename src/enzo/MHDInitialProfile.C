#include <string.h>
#include <stdio.h>
#include "myenzoutils.h"
#include "DebugMacros.h"
#include "DebugTools.h"
#include "MHDInitialProfile.h"
#include "ErrorExceptions.h"

#define MYSWHITESPACE(c) (c==' ' || c=='\t' || c=='\r' || c=='\n')

/**
 * Returns the position of the next whitespace character in s following start.
 */
char* MHDInitialProfile::strSkipSpace(char* s)
{
	for(; *s != '\0'; s++)
		if(!MYSWHITESPACE(*s))
			break;
	return s;
}
/**
 * Returns the position of the next non-whitespace character following start.
 */
char* MHDInitialProfile::strSkipNonSpace(char* s)
{
	for(; *s != '\0'; s++)
		if(MYSWHITESPACE(*s))
			break;
	return s;
}

/**
 * Returns the count of non-empty words in a string separated by whitespace.
 */
long long MHDInitialProfile::countWords(char* s)
{
	long long n = 0;
	while(*(s = strSkipSpace(s)) != '\0')
	{
		n++;
		s = strSkipNonSpace(s);
	}
	return n;
}

/**
 * Returns the index of a string in an array of strings with n elements
 * or -1 if not found. NULLs in a are ignored.
 */
long long MHDInitialProfile::indexOfStr(char* s, char** a, long long n)
{
	if(s)
		for(long long i = 0; i < n; i++)
		{
			if(a[i] && !strcmp(s, a[i]))
				return i;
		}
	return -1;
}

/**
 * Returns the index of the first occurrence of an integer
 * in an array of integers or -1 if not found.
 */
long long MHDInitialProfile::indexOfInt(long long i, long long* a, long long n)
{
	for(long j = 0; j < n; j++)
		if(i == a[j])
			return j;
	return -1;
}

/**
 * Returns the would-be sorting order if nextElement was appended
 * to the array prevElements with nPrevElements number of elements,
 * which are already in prevSortingOrder.
 */
long long MHDInitialProfile::nextSortingOrder(double nextElement, double* prevElements, long long nPrevElements,
	long long prevSortingOrder)
{
	// With no previous elements this will
	// be an array with one element, nextElement.
	// We say that an array with one
	// element is sorted in ascending order.
	if(nPrevElements == 0)
		return PROFILE_ASC_SORT;

	// If the array was previously UNSORTED,
	// adding another element can't change that.
	if(prevSortingOrder == PROFILE_UNSORTED)
		return PROFILE_UNSORTED;

	// Find the sorting order of the last element and the nextElement.
	long long so = CMP_A_B(prevElements[nPrevElements - 1], nextElement);

	// If there was a single element before,
	// then we are going to have an array with two elements.
	// The sorting order of the two elements will be also
	// the sorting order of the new array, overwriting
	// the fiducial ascending order for a one-element arrays.
	if(nPrevElements == 1)
		return so;

	// We check if adding the nextElement will preserve the
	// prevSortOrder or make the array UNSORTED.
	return (prevSortingOrder == so) ? prevSortingOrder : PROFILE_UNSORTED;
}

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
 * repetitions, although this function
 * may find a valid index in an non-monotonic or unsorted array
 * for some search algorithms.
 */
long long MHDInitialProfile::findGTEIndex(double x, double xData[], long long nData, long long xSortingOrder)
{
	size_t n = nData - 1;
	size_t l = 0, m = n / 2, r = n;

	switch(xSortingOrder)
	{
	case PROFILE_ASC_SORT:
		while(l < m)
		{
			if(xData[m] <= x)
				l = m;
			else
				r = m;

			m = (l + r) / 2;
		}
		if(l == 0)
			return -(x < xData[0]);
		if(r == n)
			return r - (x < xData[n]);
		return l;
//		for(long long i = 0; i < nData; i++)
//		{
//			if(x <= xData[i])
//				return i;
//		}
//		return nData;
	case PROFILE_DESC_SORT:
		for(long long i = nData - 1; i >= 0; i--)
		{
			if(x >= xData[i])
				return i;
		}
		return nData;
	default:
		return -1;
	}
}

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
long long MHDInitialProfile::interpolate(double* y, double yData[], double x, double xData[], long long nData,
	long long xSortingOrder)
{
	if(yData == NULL)
		return -1;
	if(xData == NULL)
		return -2;
	if(xSortingOrder == PROFILE_UNSORTED)
		return -3;

	long long i0 = findGTEIndex(x, xData, nData, xSortingOrder);
	if(i0 < 0)
	{
		*y = yData[0];
		return 0;
	}

	long long i1 = i0 + 1;
	if(i1 == nRows)
	{
		*y = yData[i0];
		return 0;
	}

	double x0 = xData[i0];
	double x1 = xData[i1];
	double d = x1 - x0; // is d tiny?
	double c0 = (x1 - x) / d; // is c0 tiny?
	double c1 = (x - x0) / d; // is c1 tiny?
	*y = c0 * yData[i0] + c1 * yData[i1];
	return 0;
}

long long MHDInitialProfile::interpolateDensity(double* y, double x)
{
	if(0 < densityProfileMaxRadius && densityProfileMaxRadius < x)
		x = densityProfileMaxRadius;

	long long result = interpolate(y, densityData, x, radiusData, nRows, radiusSortingOrder);

	if(*y < densityProfileMinDensity)
		*y = densityProfileMinDensity;

	return result;
}

long long MHDInitialProfile::interpolateInternalEnergy(double* y, double x)
{
	int radiusSO = CMP_A_B(internalEnergyRadiusData[0], internalEnergyRadiusData[1]);
	return interpolate(y, internalEnergyData, x, internalEnergyRadiusData, nRowsInternalEnergy, radiusSO);
}

long long MHDInitialProfile::interpolateTemperature(double* y, double x)
{
	return interpolate(y, temperatureData, x, radiusData, nRows, radiusSortingOrder);
}

long long MHDInitialProfile::interpolateRadialVelocity(double* y, double x)
{
	return interpolate(y, radialVelocityData, x, radiusData, nRows, radiusSortingOrder);
}

/**
 * Returns the type of the line being processed.
 */
/**
 * Parses a string of whitespace-separated words.
 * colNames is used to store pointers to the words.
 * Each word is allocated space dynamically and terminated by '\0'.
 * Returns the count.
 */
long long MHDInitialProfile::parseColumnNames(char* line, char** colNames)
{
	long long nCols = 0;
	char* wordend = line;

	while(*(line = strSkipSpace(wordend)) != '\0')
	{
		wordend = strSkipNonSpace(line);
		long long len = wordend - line;
		colNames[nCols] = new char[len + 1];
		strncpy(colNames[nCols], line, len);
		colNames[nCols][len] = '\0';
		nCols++;
	}
	return nCols;
}

void MHDInitialProfile::free()
{
	delete[] colData;
	delete[] colNames;
	delete[] colNamesToKeep;
	delete[] colSortingOrders;
	colData          = NULL;
	colNames         = NULL;
	colNamesToKeep   = NULL;
	colSortingOrders = NULL;
}

/**
 * Profile::init()
 *
 * Initializes field members to reflect an empty profile with no data
 * nor columns. Used by constructors and destructor.
 */
void MHDInitialProfile::init()
{
	nRows = nRowsAllocated = nCols = nColsToKeep = 0;
	nRowsAllocateFirst = nRowsAllocateInc = 1024;
	colData = NULL;
	colNames = colNamesToKeep = NULL;
	colSortingOrders = colNumsToKeep = NULL;
	radiusColumnName = densityColumnName = internalEnergyColumnName = temperatureColumnName = radialVelocityColumnName =
	NULL;
	radiusData = densityData = internalEnergyData = temperatureData = radialVelocityData = NULL;
	radiusIndex = densityIndex = internalEnergyIndex = temperatureIndex = radialVelocityIndex = -1;
	radiusSortingOrder = PROFILE_UNSORTED;
	requestedTime = frameTimeFound = 0;
}

void MHDInitialProfile::addColToKeep(char* colName, char** thisColName)
{
	if(colName && *colName == '\0')
		colName = NULL;
	if(thisColName)
		*thisColName = colName;
	if(colName == NULL)
		return;

	if(colNamesToKeep)
		colNamesToKeep[nColsToKeep] = colName;

	nColsToKeep++;
}

void MHDInitialProfile::init(char* radiusColumnName, char* densityColumnName, char* internalEnergyColumnName,
	char* temperatureColumnName, char* radialVelocityColumnName, double DensityProfileMaxRadius,
	double DensityProfileMinDensity)
{
	init();

	this->densityProfileMinDensity = DensityProfileMinDensity;
	this->densityProfileMaxRadius = DensityProfileMaxRadius;

	for(int pass = 0; pass < 2; pass++)
	{
		/*
		 * addColsToKeep increases nColtoKeep and assigns thisColName as necessary.
		 * The first pass counts the column names to keep and allocates colNamesToKeep.
		 * The second pass populates the elements of colNamesToKeep.
		 */
		nColsToKeep = 0;
		addColToKeep(radiusColumnName, &this->radiusColumnName);
		addColToKeep(radialVelocityColumnName, &this->radialVelocityColumnName);
		addColToKeep(densityColumnName, &this->densityColumnName);
		addColToKeep(temperatureColumnName, &this->temperatureColumnName);
		addColToKeep(internalEnergyColumnName, &this->internalEnergyColumnName);

		if(pass == 0)
			arr_newset(&colNamesToKeep, nColsToKeep, NULL);
	}
}

/**
 * Allocates new colDolata, colNames and colSortingOrders
 * initialized with NULLs.
 * Doesn't allocate any data rows. See also: allocateMoreRows(...).
 */
long long MHDInitialProfile::allocateCols(long long nCols)
{
	this->nCols = nCols;
	arr_newset(&colData, nCols, NULL);
	arr_newset(&colNames, nCols, NULL);
	arr_newset(&colSortingOrders, nCols, PROFILE_UNSORTED);
	return nCols;
}

/**
 * Increases the allocation of the data rows arrays by a given number,
 * preserving the data that is already in the arrays.
 */
long long MHDInitialProfile::allocateMoreRows(long long nMoreRows)
{
	if(nMoreRows <= 0)
	{
		nMoreRows = (nRowsAllocated > 0) ? (nRowsAllocateInc) : (nRowsAllocateFirst);
	}
	if(nMoreRows <= 0)
		nMoreRows = 1024;
	long long nNewRows = nRowsAllocated + nMoreRows;
//	printf("Increasing number of allocated rows from %d to %d.\n", nRowsAllocated, nNewRows);
	for(long long i = 0; i < nCols; i++)
	{
		if(colNames[i] == NULL)
			continue;

		double* dest = new double[nNewRows];
		double* src = colData[i];
		if(src != NULL && nRows > 0)
		{
			memcpy((void*) dest, (const void*) src, sizeof(*src) * nRows);
			delete[] src;
		}
		colData[i] = dest;
	}

//	printf("%d rows allocated.\n", nNewRows);
	return nRowsAllocated = nNewRows;
}

/**
 * Returns the index of a column name if existing or -1.
 */
long long MHDInitialProfile::findColIndex(char* colName)
{
	return indexOfStr(colName, colNames, nCols);
}

/**
 * Returns a pointer to the column data if index and name exist.
 */
double* MHDInitialProfile::findCol(long long colIndex)
{
	return ((0 <= colIndex) && (colIndex < nCols)) ? colData[colIndex] : NULL;
}

/**
 * Returns a pointer to the column data if index and name exist.
 */
double* MHDInitialProfile::findCol(char* colName)
{
	return findCol(findColIndex(colName));
}

/**
 * Interpolates y(x) bsaed on the column indices.
 */
long long MHDInitialProfile::interpolate(double* y, long long yColIndex, double x, long long xColIndex)
{
	if(yColIndex < 0 || nCols <= yColIndex)
		return -1;
	if(xColIndex < 0 || nCols <= xColIndex)
		return -2;
	if(colSortingOrders[xColIndex] == PROFILE_UNSORTED)
		return -3;
	return interpolate(y, colData[yColIndex], x, colData[xColIndex], nRows, colSortingOrders[xColIndex]);
}

/**
 * Interpolates y(x) based on column names.
 */
long long MHDInitialProfile::interpolate(double* y, char* yname, double x, char* xname)
{
	return interpolate(y, findColIndex(yname), x, findColIndex(xname));
}

/**
 * Turns column names in colNames into NULL,
 * if not present in the colNamesToKepp or
 * the colNumsToKeep lists.
 */
long long MHDInitialProfile::removeCols()
{
	if(colNamesToKeep != NULL)
	{
		for(long long i = 0; i < nCols; i++)
			if(0 > indexOfStr(colNames[i], colNamesToKeep, nColsToKeep))
			{
//				printf("Column '%s' found, ignored.\n", p->colNames[i]);
				colNames[i] = NULL;
			}
			else
			{
//				printf("Column '%s' kept.\n", p->colNames[i]);
			}
	}
	else if(colNumsToKeep != NULL)
	{
		for(long long i = 0; i < nCols; i++)
			if(0 > indexOfInt(i, colNumsToKeep, nColsToKeep))
				colNames[i] = NULL;
	}

	long long nRemaining = 0;
	for(long long i = 0; i < nCols; i++)
		nRemaining += colNames[i] != NULL;
	return nRemaining;
}

void MHDInitialProfile::identifyNamedCols()
{
	if(0 <= (radiusIndex = findColIndex(radiusColumnName)))
	{
		radiusData = colData[radiusIndex];
		radiusSortingOrder = colSortingOrders[radiusIndex];
	}
	if(0 <= (radialVelocityIndex = findColIndex(radialVelocityColumnName)))
		radialVelocityData = colData[radialVelocityIndex];
	if(0 <= (densityIndex = findColIndex(densityColumnName)))
		densityData = colData[densityIndex];
	if(0 <= (temperatureIndex = findColIndex(temperatureColumnName)))
		temperatureData = colData[temperatureIndex];
	if(0 <= (internalEnergyIndex = findColIndex(internalEnergyColumnName)))
		internalEnergyData = colData[internalEnergyIndex];
	if(internalEnergyData)
		internalEnergyRadiusData = radiusData;
}

/**
 * Parses a columns header line to allocate the correct number of columns
 * and save the column names. Returns the number of columns retained,
 * which should be equal to nColsToKeep.
 */
long long MHDInitialProfile::processColNamesLine(char* line)
{
	long long n = 0;
	if(0 >= (n = countWords(line)))
	{
		printf("Error: %d words counted in '%s'.\n", n, line);
		return n;
	}
	if(0 >= (n = allocateCols(n)))
	{
		printf("Error: %d columns allocated in '%s'.\n", n, line);
		return n;
	}
	if(0 >= (n = parseColumnNames(line, colNames)))
	{
		printf("Error: %d column names parsed in '%s'.\n", n, line);
		return n;
	}
	long long m;
	arr_printf_csv("Column names found: ", "%s", colNames, nCols);
	if(0 >= (m = removeCols()))
	{
		printf("Error: %d data columns retained to be read. '%s'\n", m, line);
	}
	arr_printf_csv("Column names retained:", "%s", colNames, nCols);
	return nCols = n;
}

/**
 * Parses a data line as a whitespace-separated list of floats.
 * Stores the data in the data rows arrays and increments nDataRows.
 * Increases the allocation for data rows if necessary.
 * Columns with NULL names are skipped.
 */
long long MHDInitialProfile::processDataLine(char* line, long long* nLinesWithExtraCols)
{
	long long nWords = 0;
	while(*(line = strSkipSpace(line)) != '\0')
	{
		if(nWords >= nCols)
		{
			*nLinesWithExtraCols++;
			break;
		}

		if(colNames[nWords] != NULL)
		{
			double x;
			long long nScanned = sscanf(line, "%lf", &x);
			if(nScanned != 1)
			{
				printf("Problem parsing float '%s'", line);
				continue;
			}

			if(nRows + 1 > nRowsAllocated)
				allocateMoreRows(0);
			double* data = colData[nWords];
			data[nRows] = x;
			colSortingOrders[nWords] = nextSortingOrder(x, data, nRows, colSortingOrders[nWords]);
		}
		nWords++;
		line = strSkipNonSpace(line);
	}
	nRows++;
	return nWords;
}

/**
 * Sets to NULL column names in varNames if they are not in the colNamesToKeep
 * list. If colNamesToKeep is not NULL the colNumsToKeep is ignored. Otherwise
 * sets to NULL col names
 * @param colNamesToKeep
 * @param colNumsToKeep
 * @param nColsToKeep
 * @return
 */
long long MHDInitialProfile::lineType(char* line)
{
	// Skip over leading spaces and tabs.
	long long a = 0;
	line = strSkipSpace(line);

	switch(*line)
	{
	case '\0':
		return PROFILE_EMPTY_LINE; // Reached the end of the string.
	case '#':
		return PROFILE_COMMENT_LINE; // The first non-whitespace is a hash.
		// If the line begins with a number, i.e. starts with a
		// plus, '+', minus, '-', a period, '.', or a digit ('0'..'9')
		// treat it as a data row.
	case '+':
	case '-':
	case '.':
	case '0':
	case '1':
	case '2':
	case '3':
	case '4':
	case '5':
	case '6':
	case '7':
	case '8':
	case '9':
		return PROFILE_DATA_LINE;

		// Any other is a columns header:
	default:
		return PROFILE_COL_NAMES_LINE;
	}
}

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
long long MHDInitialProfile::lineTypePAH01(char* line)
{
	long long lt = lineType(line);
	if(lt != PROFILE_DATA_LINE)
		return lt;
	switch(countWords(line))
	{
	case 4:
		return PROFILE_TIME_LINE;
	default:
		return PROFILE_DATA_LINE;
	}
}

long long MHDInitialProfile::read(char* filename, char* format, double atTime)
{
	long long n;
	if(!strcmp(format, "PAH01"))
	{
		n = readPAH01(filename, atTime);
	}
	else if(!strcmp(format, "PAH02"))
	{
		n = readPAH02(filename, atTime);
	}
	else
	{
		char s[1024];
		sprintf(s, "Invalid profile format '%s'", format);
		throw(EnzoFatalException(s, __FILE__, __LINE__));
	}
	identifyNamedCols();
	TRACEF("Profile successfully read: %lld data rows in %lld columns.\n", nRows, nCols);
	return n;
}

/**
 * A factory method for the PAH01 format. Creates a ProfilePAH01
 * from the data in a given file, taking the profiles at the first
 * time that exceeds or is equal to the requested *time*.
 */
long long MHDInitialProfile::readPAH01(char* filename, double atTime)
{
	this->requestedTime = atTime;
	FILE *file;
	TRACEF("Opening profile'%s'\n", filename);
	if((file = fopen(filename, "r")) == NULL)
		return -1;

	char colNames[] = "mass dmass radius density temperature gamma col7 col8";
	processColNamesLine(colNames);

	for(long long i = 0; i < nCols; i++)
	{
		if(colNames[i] == NULL)
			printf("Col %d will be ignored.\n", i);
		else
			printf("Col %d, \"%s\" will be parsed\n", i, colNames[i]);
	}

	long long lineNum = 0;
	long long nLinesWithExtraCols = 0;
	bool searchingTime = true;
	long long nColHeadersProcessed = 0;

	const size_t MAX_LINE_LENGTH = 4096;
	char line[MAX_LINE_LENGTH];
	while(fgets(line, MAX_LINE_LENGTH, file))
	{
		lineNum++;
		if(strlen(line) == MAX_LINE_LENGTH - 1)
		{
			char s[256];
			snprintf(s, 256, "Line %lld in '%s' possibly exceeds the maximum line length, %lld. Aborting...", lineNum,
						filename, MAX_LINE_LENGTH);
			EnzoFatalException(s, "MHDInitialProfile::readPAH01", 0);
		}
		switch(lineTypePAH01(line))
		// Use 'continue' to read the next line.
		// Use 'break' to exit the switch and quit reading.
		{
		case PROFILE_EMPTY_LINE:
		case PROFILE_COMMENT_LINE:
			continue;
		case PROFILE_COL_NAMES_LINE:
//			if (searchingTime) continue;
//			if (nColHeadersProcessed) break;
//			if (0 < profileProcessColNamesLine(line, p)) nColHeadersProcessed++;
			continue;
		case PROFILE_TIME_LINE:
		{
			if(!searchingTime)
				break; //This means we already found a time header with matching time
			// and most likely read all the data below it and now we reached the next time header.

			long long nRows;
			double t;
			long long nScanned = sscanf(line, "%lld %lf", &nRows, &t);
			if(nScanned != 2)
			{
				printf("problem parsing time header %s\n", line);
			}
			searchingTime = (requestedTime > t);
			if(!searchingTime)
			{
				frameTimeFound = requestedTime;
				TRACEF("Time header for %f found on line num %d. %d rows\n", requestedTime, lineNum, nRows);
				nRowsAllocateFirst = nRows;
			}
			continue;
		}
		case PROFILE_DATA_LINE:
			if(searchingTime)
				continue;
//			if (1 > nColHeadersProcessed) continue;
//			if (1 < nColHeadersProcessed) break;
			processDataLine(line, &nLinesWithExtraCols);
			if(nLinesWithExtraCols > 0)
			{
				printf("Warning: More data columns than column names found on %d lines.\n", nLinesWithExtraCols);
			}
			continue;
		default:
			continue;
		}

		break;
	}

	fclose(file);
	printf("Reading finished.\n");
	if(searchingTime)
	{
		printf("Time header for %f not found.\n", requestedTime);
		return -1;
	}

//	printf("%d data lines in %d cols.\n", p->nRows, p->nCols);
//	for (long long j = 0; j < p->nCols; j++)
//	{
//		printf("%p ", p->colData[j]);
//	}
//	printf("\n");
//	for (long long j = 0; j < p->nCols; j++)
//	{
//		printf("%p ", p->colNames[j]);
//	}
//	printf("\n");
//	for (long long i = 0; i < 10; i++)
//	{
//		printf("%3d: ", i);
//		for (long long j = 0; j < p->nCols; j++)
//		{
//			if (p->colNames[j] != NULL)
//			{
//				printf("(%d,%d)%f ", i,j,p->colData[j][i]);
//			}
//			else
//			{
//				printf("[%d,%d]%f ", i,j,0.0);
//			}
//		}
//		printf("\n");
//	}
//	ENZO_FAIL("after profile dump");
	return 0;
}

/**
 * Overridden to distinguish between a TIME_LINE and a PROFILE_LINE,
 * both of which look like DATA_LINEs to the base method.
 */
long long MHDInitialProfile::lineTypePAH02(char* line)
{
	long long lt = lineType(line);
	if(lt != PROFILE_DATA_LINE)
		return lt;
	switch(countWords(line))
	{
	case 3:
		return PROFILE_TIME_LINE;
	default:
		return PROFILE_DATA_LINE;
	}
}

long long MHDInitialProfile::readPAH02(char* filename, double atTime)
{
	this->requestedTime = atTime;

	FILE *file;
	TRACEF("Opening profile PAH02 '%s'", filename);
	if((file = fopen(filename, "r")) == NULL)
	{
		printf("Can't open '%s'.\n", filename);
		return -1;
	}
//	printf("PAH02 open\n");

//	for (long long i = 0; i < p->nCols; i++)
//	{
//		if (p->colNames[i] == NULL) printf("Col %d will be ignored.\n", i);
//		else printf("Col %d, \"%s\" will be parsed\n", i, p->colNames[i]);
//	}

	long long lineNum = 0;
	long long nLinesWithExtraCols = 0;
	bool searchingTime = true;
	long long nColHeadersProcessed = 0;

	const size_t MAX_LINE_LENGTH = 4096;
	char line[MAX_LINE_LENGTH];
	while(fgets(line, MAX_LINE_LENGTH, file))
	{
		lineNum++;
		if(strlen(line) == MAX_LINE_LENGTH - 1)
		{
			char s[256];
			snprintf(s, 256, "Line %lld in '%s' possibly exceeds the maximum line length, %lld. Aborting...", lineNum,
						filename, MAX_LINE_LENGTH);
			EnzoFatalException(s, "MHDInitialProfile::readPAH01", 0);
		}

		switch(lineTypePAH02(line))
		// Use 'continue' to read the next line.
		// Use 'break' to exit the switch and quit reading.
		{
		case PROFILE_EMPTY_LINE:
		case PROFILE_COMMENT_LINE:
			continue;
		case PROFILE_COL_NAMES_LINE:
			if(searchingTime)
				continue;
			if(nColHeadersProcessed)
				break;
			if(0 < processColNamesLine(line))
				nColHeadersProcessed++;
			continue;
		case PROFILE_TIME_LINE:
		{
			if(!searchingTime)
				break; //This means we already found a time header with matching time
			// and most likely read all the data below it and now we reached the next time header.

			double t;
			long long nScanned = sscanf(line, "%lf", &t);
			if(nScanned != 1)
			{
				printf("problem parsing time header %s\n", line);
			}
			searchingTime = (requestedTime > t);
			if(!searchingTime)
			{
				TRACEF("Time header for t=%f found on line num %d. %d rows.", requestedTime, lineNum);
			}
			continue;
		}
		case PROFILE_DATA_LINE:
			if(searchingTime)
				continue;
			if(1 > nColHeadersProcessed)
				continue;
			if(1 < nColHeadersProcessed)
				break;
			processDataLine(line, &nLinesWithExtraCols);
			if(nLinesWithExtraCols > 0)
			{
				printf("Warning: More data columns than column names found on %d lines.\n", nLinesWithExtraCols);
			}
			continue;
		default:
			continue;
		}

		break;
	}

	fclose(file);
//	printf("Reading profile finished (%s).\n", filename);
	if(searchingTime)
	{
		printf("Time header for %f not found.\n", requestedTime);
		return -1;
	}

//	printf("%d data lines in %d cols.\n", p->nRows, p->nCols);
//	for (long long j = 0; j < p->nCols; j++)
//	{
//		printf("%p ", p->colData[j]);
//	}
//	printf("\n");
//	for (long long j = 0; j < p->nCols; j++)
//	{
//		printf("%p ", p->colNames[j]);
//	}
//	printf("\n");
//	for (long long i = 0; i < 10; i++)
//	{
//		printf("%3d: ", i);
//		for (long long j = 0; j < p->nCols; j++)
//		{
//			if (p->colNames[j] != NULL)
//			{
//				printf("(%d,%d)%f ", i,j,p->colData[j][i]);
//			}
//			else
//			{
//				printf("[%d,%d]%f ", i,j,0.0);
//			}
//		}
//		printf("\n");
//	}
//	ENZO_FAIL("after profile dump");
	return 0;
}

