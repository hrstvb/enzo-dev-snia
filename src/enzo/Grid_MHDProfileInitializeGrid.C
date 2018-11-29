/***********************************************************************
 /
 /  GRID CLASS (INITIALIZE THE MHD BLAST GRID)
 /
 /  written by: David Collins
 /  date:       2004-2013
 /
 /  modified1:  Boyan Hristov
 /  date:       2015
 /  		Added perturbation method 33 creating a plane interface
 /		between zones A and B perturbed by a sine offset.
 /
 /  PURPOSE:  See MHDBlastInitialize.C for parameters
 /
 /  RETURNS:
 /    SUCCESS or FAIL
 /
 ************************************************************************/

#include <math.h>
//#include <stdio.h>
//#include <stdlib.h>
#include <string.h>
#include <cctype> //isspace

#include "myenzoutils.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "phys_constants.h"

//using namespace std;

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

int MakeFieldConservative(int field);
int FindField(int field, int farray, int numfields);
int MHDProfileInitExactB(float* Bx, float* By, float* Bz, FLOAT x, FLOAT y, FLOAT z);

/**
 * Returns the position of the next whitespace character in s following start.
 */
char* profileStrSkipSpace(char* s)
{
	for(; *s != '\0'; s++)
		if(!isspace(*s))
			break;
	return s;
}
/**
 * Returns the position of the next non-whitespace character following start.
 */
char* profileStrSkipNonSpace(char* s)
{
	for(; *s != '\0'; s++)
		if(isspace(*s))
			break;
	return s;
}

/**
 * Returns the count of non-empty words in a string separated by whitespace.
 */
int profileCountWords(char* s)
{
	int n = 0;
	while(*(s = profileStrSkipSpace(s)) != '\0')
	{
		n++;
		s = profileStrSkipNonSpace(s);
	}
	return n;
}

/**
 * Returns the index of a string in an array of strings with n elements
 * or -1 if not found. NULLs in a are ignored.
 */
int profileIndexOfStr(char* s, char** a, int n)
{
	if(s == NULL)
		return -1;
	for(int i = 0; i < n; i++)
	{
		if(a[i] == NULL)
			continue;
		if(!strcmp(s, a[i]))
			return i;
	}
	return -1;
}

/**
 * Returns the index of the first occurrence of an integer
 * in an array of integers or -1 if not found.
 */
int profileIndexOfInt(int i, int* a, int n)
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
int profileNextSortingOrder(float nextElement, float* prevElements, int nPrevElements, int prevSortingOrder)
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
	int so = CMP_A_B(prevElements[nPrevElements - 1], nextElement);

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
 * repetitions, i.e. strictly monotonic, even though this function
 * may find a valid index in an non-monotonic or unsorted array
 * for some search algorithms.
 */
int profileFindGTEIndex(float x, float xData[], int nData, int xSortingOrder)
{
	switch(xSortingOrder)
	{
	case PROFILE_ASC_SORT:
		for(int i = 0; i < nData; i++)
		{
			if(x <= xData[i])
				return i;
		}
		return nData;
	case PROFILE_DESC_SORT:
		for(int i = nData - 1; i >= 0; i--)
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
int profileInterpolate(float* y, float yData[], float x, float xData[], int nData, int xSortingOrder)
{
	int i1 = profileFindGTEIndex(x, xData, nData, xSortingOrder);
	if(i1 < 0)
		return i1;

	if(i1 == 0)
	{
		*y = yData[0];
		return 0;
	}

	int i0 = i1 - 1;
	if(i1 == nData)
	{
		*y = yData[i0];
		return 0;
	}

	float x0 = xData[i0];
	float x1 = xData[i1];
	float d = x1 - x0; // is d tiny?
	float c0 = (x1 - x) / d; // is c0 tiny?
	float c1 = (x - x0) / d; // is c1 tiny?
	*y = c0 * yData[i0] + c1 * yData[i1];
	return 0;
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
int profileParseColumnNames(char* line, char** colNames)
{
	int nCols = 0;
	char* wordend = line;

	while(*(line = profileStrSkipSpace(wordend)) != '\0')
	{
		wordend = profileStrSkipNonSpace(line);
		int len = wordend - line;
		colNames[nCols] = new char[len + 1];
		strncpy(colNames[nCols], line, len);
		nCols++;
	}

	return nCols;
}

typedef struct
{
	int nCols, nRows, nColsToKeep, nRowsAllocated, nRowsAllocateFirst, nRowsAllocateInc;
	char **colNames, **colNamesToKeep;int *colSortingOrders, *colNumsToKeep;float time;float **colData;
} profilestruct;

/**
 * The destructor.
 */
void profileFree(profilestruct* p)
{
	delete[] p->colData;
	delete[] p->colNames;
	delete[] p->colSortingOrders;
}

/**
 * Profile::init()
 *
 * Initializes field members to reflect an empty profile state with no data
 * nor columns. Used by constructors and destructor.
 */
void profileInit(profilestruct* p)
{
	p->nRows = p->nRowsAllocated = p->nCols = p->nColsToKeep = 0;
	p->nRowsAllocateFirst = p->nRowsAllocateInc = 1024;
	p->colData = NULL;
	p->colNames = p->colNamesToKeep = NULL;
	p->colSortingOrders = p->colNumsToKeep = NULL;
	p->time = 0;
}

/**
 * Allocates new colDolata, colNames and colSortingOrders
 * initialized with NULLs.
 * Doesn't allocate any data rows. See also: allocateMoreRows(...).
 */
int profileAllocateCols(int nCols, profilestruct* p)
{
	p->nCols = nCols;
	p->colData = new float*[nCols];
	p->colNames = new char*[nCols];
	p->colSortingOrders = new int[nCols];

	for(int i = 0; i < nCols; i++)
	{
		p->colData[i] = NULL;
		p->colNames[i] = NULL;
		p->colSortingOrders[i] = 0;
	}

	return nCols;
}

/**
 * Increases the allocation of the data rows arrays by a given number,
 * preserving the data that is already in the arrays.
 */
int profileAllocateMoreRows(int nMoreRows, profilestruct* p)
{
	if(nMoreRows <= 0)
	{
		nMoreRows = (p->nRowsAllocated > 0) ? (p->nRowsAllocateInc) : (p->nRowsAllocateFirst);
	}
	if(nMoreRows <= 0)
		nMoreRows = 1024;
	int nNewRows = p->nRowsAllocated + nMoreRows;
	printf("Increasing number of allocated rows from %d to %d.\n", p->nRowsAllocated, nNewRows);
	for(int i = 0; i < p->nCols; i++)
	{
		if(p->colNames[i] == NULL)
			continue;

		float* dest = new float[nNewRows];
		float* src = p->colData[i];
		if(src != NULL && p->nRows > 0)
		{
			memcpy((void*) dest, (const void*) src, sizeof(*src) * p->nRows);
			delete[] src;
		}
		p->colData[i] = dest;
	}

	printf("%d rows allocated.\n", nNewRows);
	return p->nRowsAllocated = nNewRows;
}

/**
 * Returns the index of a column name if existing or -1.
 */
int profileFindColIndex(char* colName, profilestruct* p)
{
	return profileIndexOfStr(colName, p->colNames, p->nCols);
}

/**
 * Interpolates y(x) bsaed on the column indices.
 */
int profileInterpolate(float* y, int yColIndex, float x, int xColIndex, profilestruct* p)
{
	if(yColIndex < 0 || p->nCols <= yColIndex)
		return -1;
	if(xColIndex < 0 || p->nCols <= xColIndex)
		return -2;
	if(p->colSortingOrders[xColIndex] == PROFILE_UNSORTED)
		return -3;
	return profileInterpolate(y, p->colData[yColIndex], x, p->colData[xColIndex], p->nRows,
								p->colSortingOrders[xColIndex]);
}

/**
 * Interpolates y(x) based on column names.
 */
int profileInterpolate(float* y, char* yname, float x, char* xname, profilestruct* p)
{
	return profileInterpolate(y, profileFindColIndex(yname, p), x, profileFindColIndex(xname, p), p);
}

/**
 * Turns column names in colNames into NULL,
 * if not present in the colNamesToKepp or
 * the colNumsToKeep lists.
 */
int profileRemoveCols(profilestruct* p)
{
	if(p->colNamesToKeep != NULL)
	{
		for(int i = 0; i < p->nCols; i++)
			if(0 > profileIndexOfStr(p->colNames[i], p->colNamesToKeep, p->nColsToKeep))
				p->colNames[i] = NULL;
	}
	else if(p->colNumsToKeep != NULL)
	{
		for(int i = 0; i < p->nCols; i++)
			if(0 > profileIndexOfInt(i, p->colNumsToKeep, p->nColsToKeep))
				p->colNames[i] = NULL;
	}

	int nRemaining = 0;
	for(int i = 0; i < p->nCols; i++)
		nRemaining += p->colNames[i] != NULL;

	return nRemaining;
}

/**
 * Parses a columns header line to allocate the correct number of columns
 * and save the column names. Returns the number of columns retained,
 * which should be equal to nColsToKeep.
 */
int profileProcessColNamesLine(char* line, profilestruct* p)
{
	int n = 0;
	if(0 >= (n = profileCountWords(line)))
	{
		printf("Error: %d words counted in '%s'.\n", n, line);
		return n;
	}
	if(0 >= (n = profileAllocateCols(n, p)))
	{
		printf("Error: %d columns allocated in '%s'.\n", n, line);
		return n;
	}
	if(0 >= (n = profileParseColumnNames(line, p->colNames)))
	{
		printf("Error: %d column names parsed in '%s'.\n", n, line);
		return n;
	}
	int m;
	if(0 >= (m = profileRemoveCols(p)))
	{
		printf("Error: %d data columns retained to be read. '%s'\n", m, line);
	}
	return p->nCols = n;
}

/**
 * Parses a data line as a whitespace-separated list of floats.
 * Stores the data in the data rows arrays and increments nDataRows.
 * Increases the allocation for data rows if necessary.
 * Columns with NULL names are skipped.
 */
int profileProcessDataLine(char* line, profilestruct* p, int* nLinesWithExtraCols)
{
	int nWords = 0;
	while(*(line = profileStrSkipSpace(line)) != '\0')
	{
		if(nWords >= p->nCols)
		{
			*nLinesWithExtraCols++;
			break;
		}

		if(p->colNames[nWords] != NULL)
		{
			float x;
			int nScanned = sscanf(line, "%lf", &x);
			if(nScanned != 1)
			{
				printf("Problem parsing line '%s'", line);
				continue;
			}

			if(p->nRows + 1 > p->nRowsAllocated)
				profileAllocateMoreRows(0, p);
			float* data = p->colData[nWords];
			data[p->nRows] = x;
			p->colSortingOrders[nWords] = profileNextSortingOrder(x, data, p->nRows, p->colSortingOrders[nWords]);
		}
		nWords++;
		line = profileStrSkipNonSpace(line);
	}
	p->nRows++;
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
int profileLineType(char* line)
{
	// Skip over leading spaces and tabs.
	int a = 0;
	line = profileStrSkipSpace(line);

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
int profileLineTypePAH01(char* line)
{
	int lt = profileLineType(line);
	if(lt != PROFILE_DATA_LINE)
		return lt;
	switch(profileCountWords(line))
	{
	case 4:
		return PROFILE_TIME_LINE;
	default:
		return PROFILE_DATA_LINE;
	}
}

/**
 * A factory method for the PAH01 format. Creates a ProfilePAH01
 * from the data in a given file, taking the profiles at the first
 * time that exceeds or is equal to the requested *time*.
 */
int profileReadPAH01(char* filename, profilestruct* p)
{
	FILE *file;
	printf("Opening '%s'\n", filename);
	if((file = fopen(filename, "r")) == NULL)
		return -1;

	char colNames[] = "mass dmass radius density temperature gamma col7 col8";
	profileProcessColNamesLine(colNames, p);

	for(int i = 0; i < p->nCols; i++)
	{
		if(p->colNames[i] == NULL)
			printf("Col %d will be ignored.\n", i);
		else
			printf("Col %d, \"%s\" will be parsed\n", i, p->colNames[i]);
	}

	int lineNum = 0;
	int nLinesWithExtraCols = 0;
	bool searchingTime = true;
	int nColHeadersProcessed = 0;

	char line[MAX_LINE_LENGTH];
	while(fgets(line, MAX_LINE_LENGTH, file))
	{
		lineNum++;

		switch(profileLineTypePAH01(line))
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

			int nRows;
			float t;
			int nScanned = sscanf(line, "%lld %lf", &nRows, &t);
			if(nScanned != 2)
			{
				printf("problem parsing time header %s\n", line);
			}
			searchingTime = (p->time > t);
			if(!searchingTime)
			{
				printf("Time header for %f found on line num %d. %d rows\n", p->time, lineNum, nRows);
				p->nRowsAllocateFirst = nRows;
			}
			continue;
		}
		case PROFILE_DATA_LINE:
			if(searchingTime)
				continue;
//			if (1 > nColHeadersProcessed) continue;
//			if (1 < nColHeadersProcessed) break;
			profileProcessDataLine(line, p, &nLinesWithExtraCols);
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
		printf("Time header for %f not found.\n", p->time);
		return -1;
	}

//	printf("%d data lines in %d cols.\n", p->nRows, p->nCols);
//	for (int j = 0; j < p->nCols; j++)
//	{
//		printf("%p ", p->colData[j]);
//	}
//	printf("\n");
//	for (int j = 0; j < p->nCols; j++)
//	{
//		printf("%p ", p->colNames[j]);
//	}
//	printf("\n");
//	for (int i = 0; i < 10; i++)
//	{
//		printf("%3d: ", i);
//		for (int j = 0; j < p->nCols; j++)
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
int profileLineTypePAH02(char* line)
{
	int lt = profileLineType(line);
	if(lt != PROFILE_DATA_LINE)
		return lt;
	switch(profileCountWords(line))
	{
	case 3:
		return PROFILE_TIME_LINE;
	default:
		return PROFILE_DATA_LINE;
	}
}

int profileReadPAH02(char* filename, profilestruct* p)
{
	FILE *file;
	printf("Opening '%s'\n", filename);
	if((file = fopen(filename, "r")) == NULL)
	{
		printf("cant open\n");
		return -1;
	}
	printf("PAH02 open\n");

//	for (int i = 0; i < p->nCols; i++)
//	{
//		if (p->colNames[i] == NULL) printf("Col %d will be ignored.\n", i);
//		else printf("Col %d, \"%s\" will be parsed\n", i, p->colNames[i]);
//	}

	int lineNum = 0;
	int nLinesWithExtraCols = 0;
	bool searchingTime = true;
	int nColHeadersProcessed = 0;

	char line[MAX_LINE_LENGTH];
	while(fgets(line, MAX_LINE_LENGTH, file))
	{
		lineNum++;

		switch(profileLineTypePAH02(line))
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
			if(0 < profileProcessColNamesLine(line, p))
				nColHeadersProcessed++;
			continue;
		case PROFILE_TIME_LINE:
		{
			if(!searchingTime)
				break; //This means we already found a time header with matching time
			// and most likely read all the data below it and now we reached the next time header.

			float t;
			int nScanned = sscanf(line, "%lf", &t);
			if(nScanned != 1)
			{
				printf("problem parsing time header %s\n", line);
			}
			searchingTime = (p->time > t);
			if(!searchingTime)
			{
				printf("Time header for t=%f found on line num %d. %d rows\n", p->time, lineNum);
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
			profileProcessDataLine(line, p, &nLinesWithExtraCols);
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
	printf("Reading profile finished.\n");
	if(searchingTime)
	{
		printf("Time header for %f not found.\n", p->time);
		return -1;
	}

//	printf("%d data lines in %d cols.\n", p->nRows, p->nCols);
//	for (int j = 0; j < p->nCols; j++)
//	{
//		printf("%p ", p->colData[j]);
//	}
//	printf("\n");
//	for (int j = 0; j < p->nCols; j++)
//	{
//		printf("%p ", p->colNames[j]);
//	}
//	printf("\n");
//	for (int i = 0; i < 10; i++)
//	{
//		printf("%3d: ", i);
//		for (int j = 0; j < p->nCols; j++)
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

int grid::MHDProfileInitializeGrid(char* profileFileName, char* profileFormat, char* profileType,
									char* radiusColumnName, char* densityColumnName, char* temperatureColumnName,
									float burningTemperature,
									float burnedRadius, float profileAtTime)
{
	if(GridRank != 3)
		ENZO_FAIL("MHDProfileInitializeGrid is implemented for 3D only.")

	int useGE = DualEnergyFormalism && (HydroMethod != Zeus_Hydro); //[BH]
	if(CellWidth[0][0] <= 0)
		PrepareGridDerivedQuantities();

//  if ( PerturbMethod == 100 )
//      srand( 3449653 ); //please don't change this number.

	fprintf(stderr, "GridDim---- %d %d %d\n", GridDimension[0], GridDimension[1], GridDimension[2]);

	// Assign fieldType numbers using constants from typedefs.h,
	// as well as count the number of fields.
	NumberOfBaryonFields = 0;
	FieldType[NumberOfBaryonFields++] = Density;
	if(EquationOfState == 0)
		FieldType[NumberOfBaryonFields++] = TotalEnergy;
	if(useGE)
		FieldType[NumberOfBaryonFields++] = InternalEnergy;
	FieldType[NumberOfBaryonFields++] = Velocity1;
	FieldType[NumberOfBaryonFields++] = Velocity2;
	FieldType[NumberOfBaryonFields++] = Velocity3;
	if(WritePotential)
		FieldType[NumberOfBaryonFields++] = GravPotential;
	if(UseBurning)
		FieldType[NumberOfBaryonFields++] = Density_56Ni;   //[BH]
	if(UseMHD)
	{
		FieldType[NumberOfBaryonFields++] = Bfield1;
		FieldType[NumberOfBaryonFields++] = Bfield2;
		FieldType[NumberOfBaryonFields++] = Bfield3;
		if(HydroMethod == MHD_RK)
			FieldType[NumberOfBaryonFields++] = PhiField;
	}

	if(ProcessorNumber != MyProcessorNumber)
		return SUCCESS;

	profilestruct p;
	profileInit(&p);
	p.time = profileAtTime;
	char* colNamesToKeep[] = { radiusColumnName, densityColumnName, temperatureColumnName };
	p.colNamesToKeep = colNamesToKeep;
	p.nColsToKeep = 3;
	if(!strcmp(profileFormat, "PAH01"))
	{
		profileReadPAH01(profileFileName, &p);
	}
	else if(!strcmp(profileFormat, "PAH02"))
	{
		profileReadPAH02(profileFileName, &p);
	}
	else
	{
		char s[1024];
		sprintf(s, "Invalid profile format '%s'", profileFormat);
		ENZO_FAIL(s)
	}

	//Parameters

	int size = GetGridSize();
//  printf("Proc #%d: %d..%d, %d..%d, %d..%d;" //BH DEBUG
//          " %g..%g, %g..%g, %g..%g;" //BH DEBUG
//         " %g, %g, %g;" //BH DEBUG
//         " %g, %g, %g;" //BH DEBUG
//          "\n", MyProcessorNumber, //BH DEBUG
//         GridLeftEdge[0], GridRightEdge[0], //BH DEBUG
//         GridLeftEdge[1], GridRightEdge[1], //BH DEBUG
//         GridLeftEdge[2], GridRightEdge[2], //BH DEBUG
//         CELLCENTER(0, 0), CELLCENTER(1, 0), CELLCENTER(2, 0), //BH DEBUG
//         (GridLeftEdge[0]-CELLCENTER(0, 0))/CellWidth[0][0], //BH DEBUG
//         (GridLeftEdge[1]-CELLCENTER(1, 0))/CellWidth[1][0], //BH DEBUG
//         (GridLeftEdge[2]-CELLCENTER(2, 0))/CellWidth[2][0] //BH DEBUG
//         ); //BH DEBUG

	this->AllocateGrids();
	//
	// Hack for tests of random forcing.
	//
	if(RandomForcing == TRUE)
	{
		for(int dim = 0; dim < GridRank; dim++)
		{
			RandomForcingField[dim] = new float[size];
			for(int i = 0; i < size; i++)
				RandomForcingField[dim][i] = 1.0;
		}
	}

	int totENum, rhoNum, vxNum, vyNum, vzNum;
	int gasENum = -1, BxNum = -1, ByNum = -1, BzNum = -1, rhoNiNum = -1;
	float *totEField, *rhoField, *vxField, *vyField, *vzField;
	float *gasEField = NULL, *rhoNiField = NULL;
	float *BxField = NULL, *ByField = NULL, *BzField = NULL;

	if(IdentifyPhysicalQuantities(rhoNum, gasENum, vxNum, vyNum, vzNum, totENum, BxNum, ByNum, BzNum) == FAIL)
		ENZO_FAIL("MHDProfiileInitializeGrid: Error in IdentifyPhysicalQuantities for UseMHD==true.")
	if(UseMHD)
	{
		BxField = BaryonField[BxNum];
		ByField = BaryonField[ByNum];
		BzField = BaryonField[BzNum];
	}

	if(UseBurning)
	{
		int rhoNiNum = FindField(Density_56Ni, FieldType, NumberOfBaryonFields);
		if(rhoNiNum < 0)
			ENZO_FAIL("MHDProfiileInitializeGrid: Error in FindField for Density_56Ni.")
		rhoNiField = BaryonField[rhoNiNum];
	}

	printf("UseMHD=%lld, UseMHDCT=%lld, GridRank=%lld, size=%lld\n", UseMHD, UseMHDCT, GridRank, size);	//[BH]
	rhoField = BaryonField[rhoNum];
	vxField = BaryonField[vxNum];
	vyField = BaryonField[vyNum];
	vzField = BaryonField[vzNum];
	totEField = BaryonField[totENum];
	if(useGE)
		gasEField = BaryonField[gasENum];

	for(int f = 0; f < NumberOfBaryonFields; f++)
	{
		printf("MakeFieldConservative(%lld)==%lld", FieldType[f], MakeFieldConservative(FieldType[f]));
		if(DataLabel[f] != NULL)
			printf(" // %s\n", DataLabel[f]);
		else
			printf(" // DataLabel[%lld] not set.\n", f);
	}

	int debugnanflag = 0; //[BH]

	int radiusSO = p.colSortingOrders[profileFindColIndex(radiusColumnName, &p)];
	float* radiusData = p.colData[profileFindColIndex(radiusColumnName, &p)];
	float* densityData = p.colData[profileFindColIndex(densityColumnName, &p)];
	float* temperatureData = p.colData[profileFindColIndex(temperatureColumnName, &p)];
	float gammaMinusOne = Gamma - 1;
	for(int k = 0; k < GridDimension[2]; k++)
	{
		FLOAT z = CELLCENTER(2, k);
		FLOAT zz = square(z - SphericalGravityCenter[2]);
		for(int j = 0; j < GridDimension[1]; j++)
		{
			FLOAT y = CELLCENTER(1, j);
			FLOAT yy_zz = square(y - SphericalGravityCenter[1]) + zz;
			for(int i = 0; i <= GridDimension[0]; i++)
			{
				int index = ELT(i, j, k);
				float vx = vxField[index] = 0;
				float vy = vyField[index] = 0;
				float vz = vzField[index] = 0;
				float Bx = BxField[index] = 0;
				float By = ByField[index] = 0;
				float Bz = BzField[index] = 0;

				if(!strcmp(profileType, "RADIAL"))
				{
					FLOAT x = CELLCENTER(0, i);
					FLOAT r = sqrt(square(x - SphericalGravityCenter[0]) + yy_zz);

					float rho, T = 0;
					int retcode = profileInterpolate(&rho, densityData, r, radiusData, p.nRows, radiusSO); //g/cm**3
					//retcode = profileInterpolate(&T, temperatureData, r, radiusData, p.nRows, radiusSO); //K
					bool isBurned = (r < burnedRadius);					// || (T > 0 && T > burningTemperature));

//					pressure = EOSPolytropicFactor * rho**Gamma;
//					gasEDensity = pressure / (Gamma - 1);
//					gasE = specificGasE = gasEDensity / rho = EOSPolytropicFactor * rho**(Gamma-1) / (Gamma-1);
					float gasE = EOSPolytropicFactor * POW(rho, gammaMinusOne) / gammaMinusOne;
					if(DualEnergyFormalism)
						gasEField[index] = gasE;
					MHDProfileInitExactB(&Bx, &By, &Bz, x, y, z);
					rhoField[index] = rho;
					totEField[index] = gasE + 0.5 * (vx * vx + vy * vy + vz * vz)
							+ 0.5 * (Bx * Bx + By * By + Bz * Bz) / rho;

					if(UseBurning)
						rhoNiField[index] = (isBurned) ? rho : 0;
//					if (debug1 && j == (GridDimension[1] / 2) && k == (GridDimension[2] / 2)
//							&& i <= GridDimension[0] / 2)
//						printf("i,j,k=%03d,%04d,%04d, r=%4f, (%4f,%4f,%4f), rho=%1.3f, T=%1f, burned=%d\n", i, j, k,
//								r * 1e-5, x * 1e-5, y * 1e-5, z * 1e-5, rho * 1e-9, T * 1e-9, isBurned);
				}
			} //end baryonfield initialize
		}
	}
	profileFree(&p);

// Boiler plate code:
//  if(DualEnergyFormalism )
//    for(index=0;index<size;index++)
//    BaryonField[ gesENum ][index] =
//       BaryonField[ totENum ][index]
//      - 0.5*(BaryonField[ vNum[0] ][index]*BaryonField[ vNum[0] ][index] +
//	     BaryonField[ vNum[1] ][index]*BaryonField[ vNum[1] ][index] +
//	     BaryonField[ vNum[2] ][index]*BaryonField[ vNum[2] ][index])
//      - 0.5*(BaryonField[BxNum][index]*BaryonField[BxNum][index] +
//             BaryonField[ByNum][index]*BaryonField[ByNum][index] +
//             BaryonField[BzNum][index]*BaryonField[BzNum][index])/BaryonField[ rhoNum ][index];

	if(debugnanflag)
		printf("nans intialized."); //[BH]
	else
		printf("Initialized.\n");

	return SUCCESS;

}
