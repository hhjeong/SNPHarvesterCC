#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <windows.h>
#include <queue>
#include "Rng.h"

using namespace std;

const int numSNPs = 2000;
const int numSubjects = 4000;

vector< vector<int> > genotypes;

vector< int > label;

const int nBitLen = (numSubjects+31)/32;
unsigned int bitCases[numSNPs][3][nBitLen];
unsigned int bitControls[numSNPs][3][nBitLen];

int ans[3];
int nCase, nControl;

void getInput() {

	
    memset( bitCases, 0, sizeof bitCases );
    memset( bitControls, 0, sizeof bitControls );

	// for( int i = 0 ; i < 3 ; ++i ) { scanf("%d", &ans[i]); ans[i]--; }

	// sort( ans,ans+3);
	string line;
	getline( cin, line );
	// getline( cin, line );
    genotypes.clear();
    genotypes.resize( numSubjects, vector< int >( numSNPs ) );

    label.clear();
    label.resize( numSubjects );

    nCase = nControl = 0;
    for( int i = 0 ; i < numSubjects ; ++i ) {
        for( int j = 0 ; j < numSNPs ; ++j ) {
            scanf("%d", &genotypes[i][j] );
        }
        scanf( "%d", &label[i] );
        if( label[i] == 1 ) nCase++; else nControl++;

        for( int j = 0 ; j < numSNPs ; ++j ) {
            if( label[i] == 1 ) {
                bitCases[ j ][ genotypes[i][j] ][ i / 32 ] |= 1 << ( i%32 );
            }
            else {
                bitControls[ j ][ genotypes[i][j] ][ i / 32 ] |= 1 << ( i%32 );
            }
        }

    }
}

inline float sqr( float x ) {
	return x * x;
}

float getCHIByThreeLoci( int locusI, int locusJ, int locusK ) {

    int freq[2][27] = {0,};

    int colSum[27] = {0};
    int rowSum[2] = { nControl, nCase };


    for( int a = 0 ; a < 3 ; ++a ) for( int b = 0 ; b < 3 ; ++b ) for( int c = 0 ; c < 3 ; ++ c ) {
        for( int i = 0 ; i < nBitLen ; ++i ) {
            int howmanyCase = __popcnt( bitCases[ locusI ][ a ][ i ] & bitCases[ locusJ ][ b ][ i ] & bitCases[ locusK ][ c ][ i ] );
            int howmanyControl = __popcnt( bitControls[ locusI ][ a ][ i ] & bitControls[ locusJ ][ b ][ i ] & bitControls[ locusK ][ c ][ i] );
            int where = a * 9 + b * 3 + c;
            freq[ 1 ][ where ] += howmanyCase; 
            freq[ 0 ][ where ] += howmanyControl; 
            colSum[ where ] += howmanyCase + howmanyControl;
        }
    }

    float ret = 0.0; 
    for( int i = 0 ; i < 2 ; ++i ) {
		for( int j = 0 ; j < 27 ; ++j ) {
			float expected = rowSum[i] * colSum[j] / float(numSubjects);
			if( expected <= 1e-9 ) continue;
			ret += sqr(freq[i][j]-expected)/expected;

		}
	}

    return ret;
}


bool delSNP[numSNPs] = {0};

void init( int loci[3] ) {

	bool here[numSNPs] = {0};

	for( int i = 0 ; i < 3 ; ++i ) {
		int r = -1;
		do {
			r = RandomI(0,numSNPs-1);
		} while( delSNP[r] || here[r] );
		loci[i] = r;
		here[r] = true;
	}
}

void process() {
	
	InitWELLRNG512a();
	int numSuccessiveRun = 20;
	int numFail = 0;
   
	int iter = 0;

	float bestVal = -1;
	int bestLoci[3] = {-1,-1,-1};

	while( numFail < numSuccessiveRun ) {
		fprintf( stderr, "%d-th random run\n", ++iter);

		int loci[3] = { -1,-1,-1 };
		init( loci );

		bool indicator[numSNPs] = {0};
		for( int i = 0 ; i < 3 ; ++i ) indicator[loci[i]] = true;

		bool updated = false;
		bool hasSwap = true;

		int inner = 0;

		float curBestVal;
		while(hasSwap) {
			hasSwap = false;


			for( int i = 0 ; i < numSNPs ; ++i ) if( !indicator[i] && !delSNP[i] ) {
				curBestVal = getCHIByThreeLoci( loci[0], loci[1], loci[2] );
				int lp = -1;

				for( int j = 0 ; j < 3 ; ++j ) {
					++inner;
					int tmp = loci[j];
					loci[j] = i;

					float curVal = getCHIByThreeLoci( loci[0], loci[1], loci[2] );

					if( curBestVal < curVal ) {
						curBestVal = curVal;
						lp = j;
					}
					loci[j] = tmp;
				}

				if( lp >= 0 ) {
					fprintf( stderr, "Updated : %d,%d,%d to ", loci[0], loci[1], loci[2]);
					hasSwap = true;
					indicator[i] = true;
					indicator[loci[lp]] = false;
					loci[lp] = i;

					fprintf( stderr, "%d,%d,%d", loci[0], loci[1], loci[2]);
					fprintf( stderr, "(%f)\n", curBestVal );

				}
			}
		}
		

		for( int i = 0 ; i < 3 ; ++i ) delSNP[loci[i]] = true;

		if( bestVal < curBestVal ) {
			bestVal = curBestVal;
			bestLoci[0] = loci[0];
			bestLoci[1] = loci[1];
			bestLoci[2] = loci[2];
			numFail = 0;
		}
		else {
			numFail++;
		}

	}
   
	sort( bestLoci, bestLoci+3 );

	fprintf(stderr,"Answer : %f(%d,%d,%d)", bestVal, bestLoci[0], bestLoci[1], bestLoci[2] );
	printf("%d\t%d\t%d\t", bestLoci[0], bestLoci[1], bestLoci[2] );

	// if( bestLoci[0] == ans[0] && bestLoci[1] == ans[1] && bestLoci[2] == ans[2] ) printf("YES"); else printf("NO");

}

int main(int argc, char *argv[]) {
	// freopen("input2.txt","r",stdin);
    DWORD start = GetTickCount();
    getInput();
    process();
    DWORD end = GetTickCount();
    printf( "\t%d\n", end-start);
}


