#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "emd.h"

float dist(feature_t *F1, feature_t *F2)
{
  float dX = F1->X - F2->X, dY = F1->Y - F2->Y, dZ = F1->Z - F2->Z;
  return sqrt(dX*dX + dY*dY + dZ*dZ); 
}

int main(int argc, char **argv)
{
    if(argc < 3)
    {
      fprintf(stderr, "Use: %s <histogram file 1> <histogram file 2>\n", argv[0]);
      exit(EXIT_FAILURE);
    }

    int sizeH1, sizeH2;
    int sizeVecH1, sizeVecH2;
    int nbBinsMaxH1, nbBinsMaxH2;
    feature_t *featuresH1, *featuresH2;
    float *weightsH1, *weightsH2;

    char const* const filename1 = argv[1];
    char const* const filename2 = argv[2];
    FILE* file1 = fopen(filename1, "r");
    if(!file1)
    {
        fprintf(stderr, "File %s couldn't be opened\n", filename1);
        exit(EXIT_FAILURE);
    }
    FILE* file2 = fopen(filename2, "r");
    if(!file2)
    {
        fprintf(stderr, "File %s couldn't be opened\n", filename2);
        exit(EXIT_FAILURE);
    }
    char line[256];
    
    if(!fgets(line, sizeof(line), file1))
    {
        fprintf(stderr, "File %s is empty\n", filename1);
        exit(EXIT_FAILURE);
    }
    sizeH1=atoi(line);
    fgets(line, sizeof(line), file1);
    sizeVecH1=atoi(line); //< first two lines are number of observations and number of bins
    fgets(line, sizeof(line), file1);
    nbBinsMaxH1=atoi(line);


    if(!fgets(line, sizeof(line), file2)) //< redo
    {
        fprintf(stderr, "File %s is empty\n", filename2);
        exit(EXIT_FAILURE);
    }
    sizeH2=atoi(line);
    fgets(line, sizeof(line), file2);
    sizeVecH2=atoi(line);
    fgets(line, sizeof(line), file2);
    nbBinsMaxH2=atoi(line);

    featuresH1 = malloc(sizeVecH1*sizeof(feature_t));
    featuresH2 = malloc(sizeVecH2*sizeof(feature_t));
    weightsH1 = malloc(sizeVecH1*sizeof(float));
    weightsH2 = malloc(sizeVecH2*sizeof(float));

    int i;
    for(i=0; fgets(line, sizeof(line), file1); ++i) 
    {
        char * token;
        token = strtok (line, " ,-"); //X
        featuresH1[i].X = ((float)atoi(token))/nbBinsMaxH1;
        token = strtok (NULL, " ,-"); //Y
        featuresH1[i].Y = ((float)atoi(token))/nbBinsMaxH1;
        token = strtok (NULL, " ,-"); //Z
        featuresH1[i].Z = ((float)atoi(token))/nbBinsMaxH1;

        token = strtok (NULL, " ,-"); //Occurences
        weightsH1[i] = (float)atoi(token)/sizeH1;
    }

    fclose(file1);

    for(i=0; fgets(line, sizeof(line), file2); ++i) 
    {
        char * token;
        token = strtok (line, " ,-"); //X
        featuresH2[i].X = ((float)atoi(token))/nbBinsMaxH2;
        token = strtok (NULL, " ,-"); //Y
        featuresH2[i].Y = ((float)atoi(token))/nbBinsMaxH2;
        token = strtok (NULL, " ,-"); //Z
        featuresH2[i].Z = ((float)atoi(token))/nbBinsMaxH2;

        token = strtok (NULL, " ,-"); //Occurences
        weightsH2[i] = (float)atoi(token)/sizeH2;
    }

    fclose(file2);

    signature_t sH1 = {sizeVecH1, featuresH1, weightsH1};
    signature_t sH2 = {sizeVecH2, featuresH2, weightsH2};
    float e;

    e = emd(&sH1, &sH2, dist, 0, 0) / sqrt(3.0);

    printf("emd %s - %s %f\n", filename1, filename2, e);

    return EXIT_SUCCESS;
}

