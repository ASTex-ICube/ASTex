// classconverterKmeans.cpp : définit le point d'entrée pour l'application console.
//
// This code creates: Wang tiles from input and noise images by inpainting


#include <istream>
#include <fstream>
#include <sstream>
#include "hv/hvPicture.h"
#include <ASTex/special_io.h>
#include <ASTex/utils.h>
#include <string>
using namespace hview;

#define LOAD_FILT  // filtered image already exists


//char *fname = "animal_print_paper_7567";
//char *fname = "confiture_resized";
//char *fname = "bricks_crop2";
//char *fname = "bricks_3141233";
//char *fname = "bread2";
//char *fname = "blueberries2";
//char *fname = "berry2";
//char *fname = "blue_rust";
//char *fname = "blue_silk_2020005";
//char *fname = "bumpy_hard_concrete_texture_9261475";
//char *fname = "coarse_hair_071009";
//char *fname = "completely_rusted";
//char *fname = "concrete_170793";
//char *fname = "cracked_paint_240037";
//char *fname = "cracked_asphalt_160796";
//char *fname = "cracked_wooden_plank_7090594";
//char *fname = "creased_fabric_7458";
//char *fname = "dotted_cracked_concrete_9290255";
//char *fname = "flaked_plaster_wall_4602";
//char *fname = "/home/guingo/Documents/Code/ASTex/LRPN/images/lichen_on_stone_5140209";
//char *fname = "lined_woolen_material_2020103";
//char *fname = "melon_skin_1011130";
//char *fname = "metal_peeling_paint_4131527";
//char *fname = "mop_bottom_3020616";
//char *fname = "mosaicstones2";
//char *fname = "natural_noise_recrop";
//char *fname = "painted_flaked_concrete_222733";
//const char *fname = "/home/guingo/Documents/Code/ASTex/LRPN/images/parched_cracked_mud_rainspots_2260562";
//char *fname = "sandy_footprints_4573";
//char *fname = "scratched_then_rusted_metal_4131471";
//char *fname = "stippled_concrete_7070209";
//char *fname = "worn_moroccan_rug_with_sand_220523";

const int NCLUSTERS = 3;
const int FFTSIZE = 8;

//const double TSCALE = 0.25;
//const double TAMPX = 10.0;
//const double TAMPY = 10.0;

const int BORDER = 2; // for wang tile generation, apply blending along cuts
const int WANGBLEND = 1;
const int NTESTS = 25; // used to generate the tiles
const int NTESTSI = 20;
const int maxpatch = 4; // the higher, the more variety in tiles but more time

const int CSIZE[NCLUSTERS] = { 32 }; // used to filter image if no one exists
const int MAX_TP_LEVELS = NCLUSTERS;
const int MAX_CL = NCLUSTERS;

int main(int argc, char* argv[])
{
    if( argc < 3 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << "<out_dir>  <input dir>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string base_dir = argv[1];
    std::string filename_source = argv[2];


    std::string delimiter = ".p";
//    std::string filename_plus_ext = basename(filename_source.c_str());
//    std::string name_file = filename_plus_ext.substr(0, filename_plus_ext.find(delimiter));

	std::string filename_plus_ext = ASTex::IO::remove_path(filename_source);
	std::string name_file = ASTex::IO::remove_ext(filename_plus_ext);

    std::string path = "/home/guingo/Images/Test_methods/ImgForSynthesis/"+name_file+"/PPM/";
	const char *path_in = path.c_str();
    std::string path_out = base_dir+"/"+name_file+"/";
	const char* path_o = path_out.c_str();
	const char *fname = name_file.c_str();

    ASTex::create_directory(path_out);

    int i,j,k;

    printf("Converting Texture: %s\n", fname);
    char buff[512];
    hvPictRGB<unsigned char> pioriginal;

    printf("loading example picture...\n");
    sprintf(buff, "%s/%s.ppm",path_in, fname);
    FILE *fd = fopen(buff, "rb");
    if (fd == 0) { perror("cannot load file:"); return 1; }
    pioriginal.loadPPM(fd, 1);
    fclose(fd);
    std::cout << "Load ok" << std::endl;

    hvPictRGB<unsigned char> distances[MAX_TP_LEVELS];
    for (k=0; k<NCLUSTERS; k++)
    {
        printf("loading cluster %d...\n",k);
        sprintf(buff,"%s/%s_nclusters_%d_fftsize_%d_mask_BIN%d.ppm",path_in, fname,NCLUSTERS,FFTSIZE,k);
        FILE *fd = fopen(buff, "rb");
        if (fd == 0) { perror("cannot load file:"); return 1; }
        distances[k].loadPPM(fd, 1);
        fclose(fd);
        std::cout << " LOAD OK" << std::endl;

    }

    hvPict<unsigned short> origclasses(pioriginal.sizeX(), pioriginal.sizeY(),0);
    for (i = 0; i<pioriginal.sizeX(); i++) for (j = 0; j<pioriginal.sizeY(); j++)
        {
            int kmax = 0;
            for (k=1; k<NCLUSTERS; k++)
            {
                if (distances[kmax].get(i,j).RED()<distances[k].get(i,j).RED())
                {
                    kmax = k;
                    origclasses.update(i, j, (unsigned short)k);
                    //if (i % 2 == 0 && j % 2 == 0 && i / 2<classes[LEVEL].sizeX() && j / 2<classes[LEVEL].sizeY()) classes[LEVEL].update(i / 2, j / 2, (unsigned short)k);
                    //classescol[LEVEL].update(i,j,hvColRGB<unsigned char>(((k*7)%127+128)%255, ((k*19)%127+128)%255, ((k*199)%127+128)%255));
                }
            }
        }
    // make class colors
    std::vector<hvColRGB<unsigned char> > lcol(MAX_CL);
    hvColRGB<double> coltab[MAX_CL];
    int npixclass[MAX_CL];
    lcol.clear();
    for (i = 0; i<NCLUSTERS; i++) { coltab[i] = hvColRGB<double>(0.0); npixclass[i] = 0; }
    for (i = 0; i<origclasses.sizeX(); i++) for (j = 0; j<origclasses.sizeY(); j++)
        {
            unsigned short ind = origclasses.get(i, j);
            if (ind >= NCLUSTERS) hvFatal("incoherent index in makeColorList");
            npixclass[ind]++;
            hvColRGB<unsigned char> col = pioriginal.get(i, j);
            coltab[ind] += (hvColRGB<double>)col;
        }
    for (i = 0; i<NCLUSTERS; i++)
        lcol.push_back(hvColRGB<unsigned char>((unsigned char)(coltab[i].RED() / (double)npixclass[i]), (unsigned char)(coltab[i].GREEN() / (double)npixclass[i]), (unsigned char)(coltab[i].BLUE() / (double)npixclass[i])));
    // save class pictures
    hvPictRGB<unsigned char> pcl(pioriginal.sizeX(), pioriginal.sizeY(), hvColRGB<unsigned char>(255));
    for (i = 0; i<pioriginal.sizeX(); i++) for (j = 0; j<pioriginal.sizeY(); j++)
        {
            k = origclasses.get(i, j);
            hvColRGB<unsigned char> cc = lcol.at(k);
            pcl.update(i, j, cc);
        }
    sprintf(buff, "%s/%s_nclusters_%d_clusters.ppm",path_o, fname, NCLUSTERS);
    fd = fopen(buff, "wb");
    if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
    pcl.savePPM(fd, 1);
    fclose(fd);

    // make filtered pictures
    hvPictRGB<unsigned char> fcl, diff, forig;
    hvPict<unsigned char> mask;

    printf("filtering ...\n");
    fcl.reset(pioriginal.sizeX(), pioriginal.sizeY(), hvColRGB<unsigned char>(255));
    forig.reset(pioriginal.sizeX(), pioriginal.sizeY(), hvColRGB<unsigned char>(255));
    diff.reset(pioriginal.sizeX(), pioriginal.sizeY(), hvColRGB<unsigned char>(255));
    mask.reset(pioriginal.sizeX(), pioriginal.sizeY(), 0);
    for (i = 0; i < pioriginal.sizeX(); i++) for (j = 0; j < pioriginal.sizeY(); j++)
        {
            k = origclasses.get(i, j);
            int ii, jj;
            for (ii = -BORDER; ii <= BORDER; ii++) for (jj = -BORDER; jj <= BORDER; jj++)
            {
                int px = i + ii, py = j + jj;
                if (px >= 0 && py >= 0 && px < pioriginal.sizeX() && py < pioriginal.sizeY())
                {
                    if (origclasses.get(px, py) != k)
                    {
                        mask.update(i, j, 255);
                    }
                }
            }
        }
        if (BORDER > 0)
        {
            mask.gaussianBlur(BORDER, BORDER); mask.gaussianBlur(BORDER, BORDER);
            mask.gaussianBlur(BORDER, BORDER); mask.gaussianBlur(BORDER, BORDER);
        }
        hvPictRGB<unsigned char> maskpi(mask, 1);
        sprintf(buff, "%s/%s_bordermask.ppm",path_o, fname);
        fd = fopen(buff, "wb");
        if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
        maskpi.savePPM(fd, 1);
        fclose(fd);

#ifdef LOAD_FILT
        //sprintf(buff, "%s_filt.ppm", fname);
        sprintf(buff,"%s/%s_Filtersize_16_FFTsize_%d_SigFreq_0.030000_STEP_5_S.ppm",path_in, fname,FFTSIZE);
        fd = fopen(buff, "rb");
        if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
        forig.loadPPM(fd, 1);
        fclose(fd);
        for (i = 0; i < pioriginal.sizeX(); i++) for (j = 0; j < pioriginal.sizeY(); j++)
        {
            hvColRGB<unsigned char> col = forig.get(i, j);
            col.blend((hvColRGB<double>)pioriginal.get(i, j), col, (double)mask.get(i, j) / 255.0);
            fcl.update(i, j, (hvColRGB<unsigned char>)col);
            hvColRGB<unsigned char> cc; cc.difference(pioriginal.get(i, j), fcl.get(i, j), 127.0, 127.0);
            diff.update(i, j, cc);
        }

#else
        for (i = 0; i<pioriginal.sizeX(); i++) for (j = 0; j<pioriginal.sizeY(); j++)
        {
            k = origclasses.get(i, j);
            int ii, jj;
            double sum = 0.0;
            hvColRGB<double> col;
            int FSIZE = CSIZE[k]/2;
            for (ii = -FSIZE; ii <= FSIZE; ii++) for (jj = -FSIZE; jj <= FSIZE; jj++)
            {
                int px = i + ii, py = j + jj;
                if (px >= 0 && py >= 0 && px < pioriginal.sizeX() && py < pioriginal.sizeY())
                {
                    if (origclasses.get(px, py) == k)
                    {
                        double weight = exp(-(double)(ii*ii + jj*jj) / (double)(FSIZE*FSIZE));
                        sum += weight;
                        hvColRGB<double> cc = (hvColRGB<double>)pioriginal.get(px, py);
                        cc *= weight;
                        col += cc;
                    }
                }
            }
            col /= sum;
            forig.update(i, j, (hvColRGB<unsigned char>)col);
            col.blend((hvColRGB<double>)pioriginal.get(i, j), col, (double)mask.get(i, j) / 255.0);
            fcl.update(i, j, (hvColRGB<unsigned char>)col);
            hvColRGB<unsigned char> cc; cc.difference(pioriginal.get(i, j), fcl.get(i, j), 127.0, 127.0);
            diff.update(i, j, cc);
        }

        sprintf(buff, "%s_filt.ppm", fname);
        fd = fopen(buff, "wb");
        if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
        forig.savePPM(fd, 1);
        fclose(fd);
#endif
        sprintf(buff, "%s/%s_filt_wb.ppm",path_o, fname);
        fd = fopen(buff, "wb");
        if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
        fcl.savePPM(fd, 1);
        fclose(fd);

        sprintf(buff, "%s/%s_diff.ppm",path_o, fname);
        fd = fopen(buff, "wb");
        if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
        diff.savePPM(fd, 1);
        fclose(fd);


        // make the noise images by inpainting
        /////////////////////////////////////////////////////////////
        hvPictRGB<unsigned char> noisei[MAX_CL], noiseisub[MAX_CL];
        hvPictRGB<unsigned char> noiseimask[MAX_CL];
        int ind;
        for (ind = 0; ind < NCLUSTERS; ind++)
        {
            noiseimask[ind].reset(pioriginal.sizeX(), pioriginal.sizeY(),hvColRGB<unsigned char>(0));
            noiseisub[ind].reset(pioriginal.sizeX(), pioriginal.sizeY(), hvColRGB<unsigned char>(255));
            hvBitmap map(pioriginal.sizeX(), pioriginal.sizeY(), false);
            for (i = 0; i < origclasses.sizeX(); i++) for (j = 0; j < origclasses.sizeY(); j++)
            {
                if (origclasses.get(i, j) == ind) map.set(i, j, true);
                if (map.get(i,j))
                {
                    if (mask.get(i,j)>0) noiseimask[ind].update(i,j,hvColRGB<unsigned char>(255-mask.get(i,j)));
                    else noiseimask[ind].update(i,j,hvColRGB<unsigned char>(255));
                    noiseisub[ind].update(i, j, diff.get(i, j));
                }
            }
            sprintf(buff, "%s/%s_nclusters_%d_noise_cl%03d.mask.inpaint.ppm",path_o, fname, NCLUSTERS, ind);
            fd = fopen(buff, "wb");
            if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
            noiseimask[ind].savePPM(fd, 1);
            fclose(fd);
            sprintf(buff, "%s/%s_nclusters_%d_noise_cl%03d.sub.ppm",path_o, fname, NCLUSTERS, ind);
            fd = fopen(buff, "wb");
            if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
            noiseisub[ind].savePPM(fd, 1);
            fclose(fd);
            hvColRGB<unsigned char> avgcol = diff.avg(map);
            map.erosion(3, 3); map.erosion(3, 3);
            if (map.count() <= 5)
                noisei[ind].reset(pioriginal.sizeX(), pioriginal.sizeY(), avgcol);
            else
            {
                noisei[ind].reset(pioriginal.sizeX(), pioriginal.sizeY(), hvColRGB<unsigned char>(127));
                printf("making noise class %d\n", ind);
                noisei[ind].lapped(diff, map, 1.0);
            }
            printf("done.\n");
            sprintf(buff, "%s/%s_nclusters_%d_noise_cl%03d.inpaint.ppm",path_o, fname, NCLUSTERS,ind);
            fd = fopen(buff, "wb");
            if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
            noisei[ind].savePPM(fd, 1);
            fclose(fd);
        }

//		exit(0);

    // make Wang tile of structure
    /////////////////////////////////////////////////////////////

        const double pow = 2.0;
        const double powfilt = 1.0;
        const int blendborder = 3;
        const int blendborderfilt = 3;
        int npatch=0;

        //double radius=(double)bradius->getValue()*0.5;
        //int maxpatch=(int)(8.0*(1.0-radius*2.0)+1.5); if (maxpatch<2) maxpatch=2;
        printf("Computing 16 Wang tiles...\n");
        //printf("maxpatch=%d, radius=%g\n", maxpatch,radius);
        //int pixrad = (int)((1.0-radius)*(double)(piexample.sizeX()<piexample.sizeY()?piexample.sizeX():piexample.sizeY()));

        hvPictRGB<unsigned char> wangpict(pioriginal.sizeX()*4, pioriginal.sizeY()*4, hvColRGB<unsigned char>(0));
        hvPictRGB<unsigned char> wangfiltpict(pioriginal.sizeX()*4, pioriginal.sizeY()*4, hvColRGB<unsigned char>(0));
        hvPictRGB<unsigned char> wangfiltwpict(pioriginal.sizeX() * 4, pioriginal.sizeY() * 4, hvColRGB<unsigned char>(0));
        hvPictRGB<unsigned char> wangnoiseimask[NCLUSTERS];
        for (ind = 0; ind < NCLUSTERS; ind++) wangnoiseimask[ind].reset(pioriginal.sizeX()*4, pioriginal.sizeY()*4, hvColRGB<unsigned char>(0));
        hvBitmap wangyet(pioriginal.sizeX()*4, pioriginal.sizeY()*4, false);
        for (i=0; i<pioriginal.sizeX()*4;i++) for (j=0; j<pioriginal.sizeY()*4; j++)
        {
            wangpict.update(i,j,pioriginal.get(i%pioriginal.sizeX(),j%pioriginal.sizeY()));
            wangfiltpict.update(i,j,forig.get(i%pioriginal.sizeX(),j%pioriginal.sizeY()));
            wangfiltwpict.update(i, j, fcl.get(i%pioriginal.sizeX(), j%pioriginal.sizeY()));
            for (ind = 0; ind < NCLUSTERS; ind++) wangnoiseimask[ind].update(i,j,noiseimask[ind].get(i%pioriginal.sizeX(),j%pioriginal.sizeY()));
        }

        hvVec2<int> fmin, fmax;
        int dimx,dimy;
        dimx = pioriginal.sizeX()/2;
        dimy = pioriginal.sizeY()/2;
        hvPict<double> imgHT, imgHB, imgVL, imgVR;
        int *pathHT = new int [dimx];
        int *pathHB = new int [dimx];
        int *pathVL = new int [dimy];
        int *pathVR = new int [dimy];

        //hvPict<hvVec3<int> > fragcenters;
        //hvPict<double> weight(pifrag.sizeX(), pifrag.sizeY(),0.0);
        //hvBitmap featborder(pifrag.sizeX(), pifrag.sizeY(), true);
        //hvBitmap feature(pifrag.sizeX(), pifrag.sizeY(), false);
        //hvBitmap tfeature(pifrag.sizeX(), pifrag.sizeY(), false);
        //hvBitmap tfragfeature(pifrag.sizeX(), pifrag.sizeY(), false);
        hvBitmap fragfeature(dimx, dimy, false);
        //pifrag.fragmentRegular(fragcenters, featborder);


        // fill from corner
        bool search=true, bordersdone=false, firstcorner=true;
        int cornercc=0;
        int tilex=0, tiley=0, tilecounter=0;
        int counter=0;
        do {
            int modulox=0; int moduloy=0;
            int locx=pioriginal.sizeX()-1,locy=pioriginal.sizeY()-1;
            search=true;

            // first manage to fill the corner and the borders
            if (!bordersdone)
            {
                const int seamlen=4;

                locx=pioriginal.sizeX()-1; locy=pioriginal.sizeY()-1;
                if (!wangyet.get(locx,locy)) search=false;
                if (search)
                {
                    locx=pioriginal.sizeX()-1; locy=pioriginal.sizeY()-1;
                    while (wangyet.get(locx,locy) && locx>pioriginal.sizeX()-4 && locy>pioriginal.sizeY()-4) { locx--; locy--; }
                    if (!wangyet.get(locx,locy)) { firstcorner=true; search=false; }
                }

                if (search) { locx=pioriginal.sizeX()/2; locy=pioriginal.sizeY()-1;
                if (!wangyet.get(locx,locy)) search=false; }
                if (search) { locx=pioriginal.sizeX()-1; locy=pioriginal.sizeY()/2;
                if (!wangyet.get(locx,locy)) search=false; }
                if (search) { locx=pioriginal.sizeX()/2; locy=2*pioriginal.sizeY()-1;
                if (!wangyet.get(locx,locy)) { search=false; moduloy=1; } }
                if (search) { locx=2*pioriginal.sizeX()-1; locy=pioriginal.sizeY()/2;
                if (!wangyet.get(locx,locy)) { search=false; modulox=1; } }

                if (search)
                {
                    locx=pioriginal.sizeX()-1; locy=pioriginal.sizeY()-1;
                    while (wangyet.get(locx,locy) && locx>1) locx--;
                    int cc=locx-1; while(cc>0 && !wangyet.get(cc,locy)) cc--;
                    if (locx>1 && (locx-cc)>seamlen) search=false;
                    else if (locx>1) for (;cc<=locx;cc++) wangyet.set(cc,locy,true);
                }

                if (search)
                {
                    locx=pioriginal.sizeX()-1; locy=pioriginal.sizeY()-1;
                    while (wangyet.get(locx,locy) && locy>1) locy--;
                    int cc=locy-1; while(!wangyet.get(locx,cc) && cc>0) cc--;
                    if (locy>1 && (locy-cc)>seamlen) search=false;
                    else if (locy>1) for (;cc<=locy;cc++) wangyet.set(locx,cc,true);
                }

                if (search)
                {
                    locx=pioriginal.sizeX()-1; locy=2*pioriginal.sizeY()-1;
                    while (wangyet.get(locx,locy) && locx>1) locx--;
                    int cc=locx-1; while(!wangyet.get(cc,locy) && cc>0) cc--;
                    if (locx>1 && (locx-cc)>seamlen) { search=false; moduloy=1; }
                    else if (locx>1) for (;cc<=locx;cc++) wangyet.set(cc,locy,true);
                }

                if (search)
                {
                    locx=2*pioriginal.sizeX()-1; locy=pioriginal.sizeY()-1;
                    while (wangyet.get(locx,locy) && locy>1) locy--;
                    int cc=locy-1; while(!wangyet.get(locx,cc) && cc>0) cc--;
                    if (locy>1 && (locy-cc)>seamlen) { search=false; modulox=1; }
                    else if (locy>1) for (;cc<=locy;cc++) wangyet.set(locx,cc,true);
                }
                while (search && locy>1)
                {
                    locx=2*pioriginal.sizeX()-1; locy=pioriginal.sizeY()-1;
                    while (wangyet.get(locx,locy) && locy>1) locy--;
                    int cc=locy-1; while(!wangyet.get(locx,cc) && cc>0) cc--;
                    if (locy>1 && (locy-cc)>seamlen) { search=false; modulox=1; }
                    else if (locy>1) for (;cc<=locy;cc++) wangyet.set(locx,cc,true);
                }

                if (search) { bordersdone=true; printf("\nborders are done...\n"); }
            }
            // manage the interior of the tiles

            if (bordersdone)
            {
                if (tilex<4 && tiley<4)
                {
                    int cc=0;
                    do {
                        cc++;
                        locx=tilex*pioriginal.sizeX()+pioriginal.sizeX()/2+int((2.0*(double)rand()/(double)RAND_MAX-1.0)*(double)(pioriginal.sizeX())/4.0);
                        locy=tiley*pioriginal.sizeY()+pioriginal.sizeY()/2+int((2.0*(double)rand()/(double)RAND_MAX-1.0)*(double)(pioriginal.sizeY())/4.0);
                        modulox=tilex; moduloy=tiley;
                    } while(cc<50 && wangyet.get(locx,locy));
                    search=false;
                }
                tilecounter++; if (tilecounter==maxpatch) { tilecounter=0; tilex++; if (tilex>=4) { tilex=0; tiley++; } }
                //search=true;
            }

            if (!search)
            {
                printf("search pos: %d,%d: ", locx,locy);
                int kk;
                hvVec2<int> minpos,minoffset;
                hvColRGB<double> minerr;
                double minf;
                for (kk=0; kk<(bordersdone?NTESTSI:NTESTS); kk++)
                {
                    double sum=0.0;
                    fmin = hvVec2<int>(int((double)pioriginal.sizeX()*(double)rand()/(double)RAND_MAX*0.5),
                            int((double)pioriginal.sizeY()*(double)rand()/(double)RAND_MAX*0.5));
                    // choose position to minimize squared difference
                    //hvVec2<int> ppos = wangsaliency.chooseWeightedMinSquareDiff(locx-3*dimx/4, locy-3*dimy/4, dimx/2, dimy/2, fmin.X(),fmin.Y(),dimx,dimy, psaliency, tfeature, weight, locx, locy, minerr);
                    hvVec2<int> ppos = wangpict.chooseMinSquareDiff(locx-dimx/2-dimx/8, locy-dimy/2-dimy/8, dimx/4, dimy/4, fmin.X(),fmin.Y(),dimx,dimy, pioriginal, minerr);

                    imgHT.squaredDifference(ppos.X(),ppos.Y()+3*dimy/4,dimx,dimy/4,wangpict,fmin.X(),fmin.Y()+3*dimy/4,pioriginal);
                    imgHT.minPathH(0,dimy/4,2,pathHT);
                    for(int ind=0; ind<dimx; ind++) sum+=imgHT.get(ind,pathHT[ind]);
                    imgHB.squaredDifference(ppos.X(),ppos.Y(),dimx,dimy/4,wangpict,fmin.X(),fmin.Y(),pioriginal);
                    imgHB.minPathH(0,dimy/4,2,pathHB);
                    for(int ind=0; ind<dimx; ind++) sum+=imgHB.get(ind,pathHB[ind]);
                    imgVR.squaredDifference(ppos.X()+3*dimx/4,ppos.Y(),dimx/4,dimy,wangpict,fmin.X()+3*dimx/4,fmin.Y(),pioriginal);
                    imgVR.minPathV(0,dimx/4,2,pathVR);
                    for(int ind=0; ind<dimy; ind++) sum+=imgVR.get(pathVR[ind],ind);
                    imgVL.squaredDifference(ppos.X(),ppos.Y(),dimx/4,dimy,wangpict,fmin.X(),fmin.Y(),pioriginal);
                    imgVL.minPathV(0,dimx/4,2,pathVL);
                    for(int ind=0; ind<dimy; ind++) sum+=imgVL.get(pathVL[ind],ind);

                    if (kk==0) { minf=sum; minpos=ppos; minoffset=fmin; }
                    else if (minf>sum) { minf=sum; minpos=ppos; minoffset=fmin;  }
                }
                hvVec2<int> startmin = minoffset;
                // refinement
                for (kk = 0; kk<NTESTSI/2; kk++)
                {
                    double sum = 0.0;
                    int fminx = startmin.X() + int((double)pioriginal.sizeX()/16.0*(1.0-2.0*(double)rand() / (double)RAND_MAX));
                    int fminy = startmin.Y() + int((double)pioriginal.sizeY() / 16.0*(1.0 - 2.0*(double)rand() / (double)RAND_MAX));
                    if (fminx < 0) fminx = 0;
                    if (fminx >= pioriginal.sizeX() / 2) fminx = pioriginal.sizeX() / 2 - 1;
                    if (fminy < 0) fminy = 0;
                    if (fminy >= pioriginal.sizeY() / 2) fminy = pioriginal.sizeY() / 2 - 1;
                    fmin = hvVec2<int>(fminx,fminy);
                    // choose position to minimize squared difference
                    //hvVec2<int> ppos = wangsaliency.chooseWeightedMinSquareDiff(locx-3*dimx/4, locy-3*dimy/4, dimx/2, dimy/2, fmin.X(),fmin.Y(),dimx,dimy, psaliency, tfeature, weight, locx, locy, minerr);
                    hvVec2<int> ppos = wangpict.chooseMinSquareDiff(locx - dimx / 2 - dimx / 8, locy - dimy / 2 - dimy / 8, dimx / 4, dimy / 4, fmin.X(), fmin.Y(), dimx, dimy, pioriginal, minerr);

                    imgHT.squaredDifference(ppos.X(), ppos.Y() + 3 * dimy / 4, dimx, dimy / 4, wangpict, fmin.X(), fmin.Y() + 3 * dimy / 4, pioriginal);
                    imgHT.minPathH(0, dimy / 4, 2, pathHT);
                    for (int ind = 0; ind<dimx; ind++) sum += imgHT.get(ind, pathHT[ind]);
                    imgHB.squaredDifference(ppos.X(), ppos.Y(), dimx, dimy / 4, wangpict, fmin.X(), fmin.Y(), pioriginal);
                    imgHB.minPathH(0, dimy / 4, 2, pathHB);
                    for (int ind = 0; ind<dimx; ind++) sum += imgHB.get(ind, pathHB[ind]);
                    imgVR.squaredDifference(ppos.X() + 3 * dimx / 4, ppos.Y(), dimx / 4, dimy, wangpict, fmin.X() + 3 * dimx / 4, fmin.Y(), pioriginal);
                    imgVR.minPathV(0, dimx / 4, 2, pathVR);
                    for (int ind = 0; ind<dimy; ind++) sum += imgVR.get(pathVR[ind], ind);
                    imgVL.squaredDifference(ppos.X(), ppos.Y(), dimx / 4, dimy, wangpict, fmin.X(), fmin.Y(), pioriginal);
                    imgVL.minPathV(0, dimx / 4, 2, pathVL);
                    for (int ind = 0; ind<dimy; ind++) sum += imgVL.get(pathVL[ind], ind);

                    if (minf>sum) { minf = sum; minpos = ppos; minoffset = fmin; }
                }
                //printf("Feature: %d,%d - %d,%d (%d,%d)\n", fmin.X(), fmin.Y(), fmax.X(), fmax.Y(),dimx, dimy);
                counter++;
                fmin=minoffset;
                imgHT.squaredDifference(minpos.X(),minpos.Y()+3*dimy/4,dimx,dimy/4,wangpict,fmin.X(),fmin.Y()+3*dimy/4,pioriginal);
                imgHT.minPathH(0,dimy/4,2,pathHT);
                imgHB.squaredDifference(minpos.X(),minpos.Y(),dimx,dimy/4,wangpict,fmin.X(),fmin.Y(),pioriginal);
                imgHB.minPathH(0,dimy/4,2,pathHB);
                imgVR.squaredDifference(minpos.X()+3*dimx/4,minpos.Y(),dimx/4,dimy,wangpict,fmin.X()+3*dimx/4,fmin.Y(),pioriginal);
                imgVR.minPathV(0,dimx/4,2,pathVR);
                imgVL.squaredDifference(minpos.X(),minpos.Y(),dimx/4,dimy,wangpict,fmin.X(),fmin.Y(),pioriginal);
                imgVL.minPathV(0,dimx/4,2,pathVL);

                // blend the patch

                bool oncorner=false;
                bool onxa=false, onxb=false, onya=false, onyb=false;
                bool pasted=false;
                fragfeature.clear(false);
                for (i=0; i<dimx; i++) for (j=0; j<dimy; j++)
                {
                    if (i>=pathVL[j] && i<=pathVR[j]+3*dimx/4 && j>=pathHB[i] && j<=pathHT[i]+3*dimy/4)
                    {
                        fragfeature.set(i,j,true);
                        int vx=minpos.X()+i;  while (vx<0) vx+=pioriginal.sizeX();  vx%=pioriginal.sizeX();
                        int vy=minpos.Y()+j;  while (vy<0) vy+=pioriginal.sizeY();  vy%=pioriginal.sizeY();
                        if ( (vx==0 || vx==pioriginal.sizeX()-1) && (vy==0 || vy==pioriginal.sizeX()-1) ) oncorner=true;
                        vx=minpos.X()+i; if (vx<0) vx+=4*pioriginal.sizeX();
                        if (vx==0 || vx==pioriginal.sizeX()-1 || vx==pioriginal.sizeX() || vx==4*pioriginal.sizeX()-1) onxa=true;
                        if (vx==2*pioriginal.sizeX()-1 || vx==2*pioriginal.sizeX() || vx==3*pioriginal.sizeX()-1 || vx==3*pioriginal.sizeX()) onxb=true;
                        vy=minpos.Y()+j; if (vy<0) vy+=4*pioriginal.sizeY();
                        if (vy==0 || vy==pioriginal.sizeY()-1 || vy==pioriginal.sizeY() || vy==4*pioriginal.sizeY()-1) onya=true;
                        if (vy==2*pioriginal.sizeY()-1 || vy==2*pioriginal.sizeY() || vy==3*pioriginal.sizeY()-1 || vy==3*pioriginal.sizeY()) onyb=true;
                    }
                }
                printf("feature %d,%d  (%d,%d) ", minpos.X(), minpos.Y(), modulox, moduloy);
                hvPict<unsigned char> pfeat(fragfeature, blendborder, 255);
                pfeat.gaussianBlur(WANGBLEND, WANGBLEND);
                hvPict<unsigned char> pfeatfilt(fragfeature, blendborderfilt, 255);
                pfeatfilt.gaussianBlur(WANGBLEND, WANGBLEND);
                if (oncorner) printf("on corner ");
                else
                {
                    if (onxa) printf("on XA "); if (onxb) printf("on XB ");
                    if (onya) printf("on YA "); if (onyb) printf("on YB ");
                }
                printf("\n");
                if (!bordersdone && (oncorner || onxa || onxb || onya || onyb) )
                {
                    minpos=hvVec2<int>(minpos.X()-modulox*pioriginal.sizeX(),minpos.Y()-moduloy*pioriginal.sizeY());
                    for (i=0; i<=4; i++) for (j=0; j<=4; j++)
                    {
                        bool doblend=false;
                        if (oncorner) { if (firstcorner) doblend=true; }
                        else
                        {
                            if (onxa && !onya && !onyb && (i==0 || i==1 || i==4)) doblend=true;
                            if (onxb && !onyb && !onya && (i==2 || i==3)) doblend=true;
                            if (onya && !onxa && !onxb && (j==0 || j==1 || j==4)) doblend=true;
                            if (onyb && !onxb && !onxa && (j==2 || j==3)) doblend=true;
                        }
                        if (doblend)
                        {
                            wangpict.blendRect(pioriginal.sizeX()*(i-1)+minpos.X(),pioriginal.sizeY()*(j-1)+minpos.Y(),fmin.X(),fmin.Y(),dimx,dimy,pioriginal,pfeat, (unsigned char)(255), pow, fragfeature,false);
                            wangfiltwpict.blendRect(pioriginal.sizeX()*(i-1)+minpos.X(),pioriginal.sizeY()*(j-1)+minpos.Y(),fmin.X(),fmin.Y(),dimx,dimy,fcl,pfeatfilt, (unsigned char)(255), powfilt, fragfeature,false);
                            wangfiltpict.blendRect(pioriginal.sizeX()*(i - 1) + minpos.X(), pioriginal.sizeY()*(j - 1) + minpos.Y(), fmin.X(), fmin.Y(), dimx, dimy, forig, pfeatfilt, (unsigned char)(255), powfilt, fragfeature, false);
                            for (ind = 0; ind < NCLUSTERS; ind++)
                                wangnoiseimask[ind].blendRect(pioriginal.sizeX()*(i-1)+minpos.X(),pioriginal.sizeY()*(j-1)+minpos.Y(),fmin.X(),fmin.Y(),dimx,dimy,noiseimask[ind],pfeatfilt, (unsigned char)(255), powfilt, fragfeature,false);
                            //wangfrag.copyWangRect(piexample.sizeX()*(i-1)+minpos.X(),piexample.sizeY()*(j-1)+minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(),pifrag,fragfeature);
                            //wangpatch.copyWangRect(piexample.sizeX()*(i-1)+minpos.X(),piexample.sizeY()*(j-1)+minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(),npatch,fragfeature);
                            //wangsaliency.copyRect(piexample.sizeX()*(i-1)+minpos.X(),piexample.sizeY()*(j-1)+minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(),psaliency,fragfeature);
                            //wangheight.blendRect(piexample.sizeX()*(i-1)+minpos.X(),piexample.sizeY()*(j-1)+minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(),heightfield,pfeathh, unsigned char(255), 1.0, ffhh);
                            wangyet.operatorOr(pioriginal.sizeX()*(i-1)+minpos.X(),pioriginal.sizeY()*(j-1)+minpos.Y(),0,0,dimx,dimy, fragfeature);
                            pasted=true;
                            cornercc=0;
                            npatch++;
                        }
                    }
                    if (oncorner && firstcorner) { firstcorner=false; cornercc=0; }
                    if (oncorner && !firstcorner) { if (cornercc>=5) firstcorner=true; else cornercc++; }
                }
                else if (bordersdone && !(oncorner || onxa || onxb || onya || onyb) )
                {

                            wangpict.blendRect(minpos.X(),minpos.Y(),fmin.X(),fmin.Y(),dimx,dimy,pioriginal,pfeat, (unsigned char)(255), pow, fragfeature,false);
                            wangfiltwpict.blendRect(minpos.X(),minpos.Y(),fmin.X(),fmin.Y(),dimx,dimy,fcl,pfeat, (unsigned char)(255), pow, fragfeature,false);
                            wangfiltpict.blendRect(minpos.X(), minpos.Y(), fmin.X(), fmin.Y(), dimx, dimy, forig, pfeat, (unsigned char)(255), pow, fragfeature, false);
                            for (ind = 0; ind < NCLUSTERS; ind++)
                                wangnoiseimask[ind].blendRect(minpos.X(),minpos.Y(),fmin.X(),fmin.Y(),dimx,dimy,noiseimask[ind],pfeat, (unsigned char)(255), pow, fragfeature,false);
                                //wangnoiseimask[ind].blendRect(pioriginal.sizeX()*(i - 1) + minpos.X(), pioriginal.sizeY()*(j - 1) + minpos.Y(), fmin.X(), fmin.Y(), dimx, dimy, noiseimask[ind], pfeat, unsigned char(255), pow, fragfeature, false);
                            //wangfrag.copyWangRect(minpos.X(),minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(),pifrag,fragfeature);
                            //wangpatch.copyWangRect(minpos.X(),minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(),npatch,fragfeature);
                            //wangsaliency.copyRect(minpos.X(),minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(),psaliency,fragfeature);
                            //wangheight.blendRect(minpos.X(),minpos.Y(),fmin.X(),fmin.Y(),fmax.X()-fmin.X(),fmax.Y()-fmin.Y(),heightfield,pfeathh, unsigned char(255), 1.0, ffhh);
                            wangyet.operatorOr(minpos.X(),minpos.Y(),0,0,dimx,dimy, fragfeature);
                            pasted=true;
                            npatch++;
                }
                else printf("feature eliminated.\n");
                if (pasted)
                {
                    sprintf(buff, "%s/%s_wang.ppm",path_o, fname);
                    fd = fopen(buff, "wb");
                    if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
                    wangpict.savePPM(fd, 1);
                    fclose(fd);
                    sprintf(buff, "%s/%s_filtw_wang.ppm",path_o, fname);
                    fd = fopen(buff, "wb");
                    if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
                    wangfiltwpict.savePPM(fd, 1);
                    fclose(fd);
                    sprintf(buff, "%s/%s_filt_wang.ppm",path_o, fname);
                    fd = fopen(buff, "wb");
                    if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
                    wangfiltpict.savePPM(fd, 1);
                    fclose(fd);
                    for (ind = 0; ind < NCLUSTERS; ind++)
                    {
                        sprintf(buff, "%s/%s_noise_cl%03d_wang.ppm",path_o, fname,ind);
                        fd = fopen(buff, "wb");
                        if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
                        wangnoiseimask[ind].savePPM(fd, 1);
                        fclose(fd);
                    }
                    //exit(0);
                }
                else if (bordersdone) tilecounter--;
            }
        } while (!search);

// Make turbulence
        /*
        sprintf(buff, "%s_filt", fname);
        hvProcWPhase proc(buff, 1, true, true, true);
        proc.doTurbulence();
        hvPictRGB<unsigned char> turb(2 * pioriginal.sizeX(), 2 * pioriginal.sizeY(), hvColRGB<unsigned char>(0));
        for (i = 0; i < turb.sizeX(); i++) for (j = 0; j < turb.sizeY(); j++)
        {
            double ppx = (double)i, ppy = (double)j;
            double tt = proc.turbulence(TSCALE*ppx, TSCALE*ppy);
            // Perlin turbulence
            //double tt = hvNoise::turbulence(0.1*(double)i, 0.1*(double)j, 2.9, 0.0001)*4.0-1.0;
            if (tt < -1.0) tt = -1.0; else if (tt > 1.0) tt = 1.0;
            turb.update(i, j, hvColRGB<unsigned char>((unsigned char)(255.0 / 2.0*(tt + 1.0))));
        }
        sprintf(buff, "%s_turb.ppm", fname);
        fd = fopen(buff, "wb");
        if (fd == 0) { perror("cannot load file:"); return 1; }
        turb.savePPM(fd, 1);
        fclose(fd);
        */

    // make the synthesis images
        /*
    hvPictRGB<unsigned char> synthi(2*pioriginal.sizeX(), 2*pioriginal.sizeY(), hvColRGB<unsigned char>());
    hvPictRGB<unsigned char> synthib(2*pioriginal.sizeX(), 2*pioriginal.sizeY(), hvColRGB<unsigned char>());
    hvPictRGB<unsigned char> synthiw(2 * pioriginal.sizeX(), 2 * pioriginal.sizeY(), hvColRGB<unsigned char>());
    int ii, jj;
    for (ii = 0; ii < 2*origclasses.sizeX(); ii++) for (jj = 0; jj < 2*origclasses.sizeY(); jj++)
    {
        //double tt = proc.turbulence(TSCALE*(double)ii, TSCALE*(double)jj);
        double ttx = proc.turbulence(TSCALE*(double)ii, TSCALE*(double)jj);
        double tty = proc.turbulence(TSCALE*(double)ii + 4.543, TSCALE*(double)jj + 9.654);
        i = ii+(int)(ttx*TAMPX); j = jj + (int)(tty*TAMPY);
        while (i < 0) i += wangfiltpict.sizeX(); while (i >= wangfiltpict.sizeX()) i -= wangfiltpict.sizeX();
        while (j < 0) j += wangfiltpict.sizeY(); while (j >= wangfiltpict.sizeY()) j -= wangfiltpict.sizeY();

        synthiw.update(ii, jj, wangpict.get(i%wangpict.sizeX(), j%wangpict.sizeY()));

        //ind = wangfiltpict.get(i, j);
        //hvColRGB<double> cc = (hvColRGB<double>)fcl.get(i, j);
        //cc.blend((hvColRGB<double>)pioriginal.get(i, j), cc, (double)mask.get(i, j) / 255.0);
        //hvColRGB<double> col = (hvColRGB<double>)noisei[ind].get(i, j);
        //col -= hvColRGB<double>(127.0);
        //col *= 1.0 - (double)mask.get(i, j) / 255.0;
        //cc += col;
        //cc.clamp(0.0, 255.0);
        //synthi.update(i, j, (hvColRGB<unsigned char>)cc);
        hvColRGB<double> cc = (hvColRGB<double>)wangfiltwpict.get(i, j);
        for (ind = 0; ind < NCLUSTERS; ind++)
        {
            int nx = ii%noisei[ind].sizeX();
            int ny = jj%noisei[ind].sizeY();
            hvColRGB<double> col = (hvColRGB<double>)noisei[ind].get(nx, ny);
            col -= hvColRGB<double>(127.0);
            col.scale((double)wangnoiseimask[ind].get(i,j).RED()/255.0);
            cc += col;
        }
        cc.clamp(0.0, 255.0);
        synthib.update(ii, jj, (hvColRGB<unsigned char>)cc);

        i = ii ; j = jj ;
        while (i < 0) i += wangfiltpict.sizeX(); while (i >= wangfiltpict.sizeX()) i -= wangfiltpict.sizeX();
        while (j < 0) j += wangfiltpict.sizeY(); while (j >= wangfiltpict.sizeY()) j -= wangfiltpict.sizeY();
        cc = (hvColRGB<double>)wangfiltwpict.get(i, j);
        for (ind = 0; ind < NCLUSTERS; ind++)
        {
            int nx = ii%noisei[ind].sizeX();
            int ny = jj%noisei[ind].sizeY();
            hvColRGB<double> col = (hvColRGB<double>)noisei[ind].get(nx, ny);
            col -= hvColRGB<double>(127.0);
            col.scale((double)wangnoiseimask[ind].get(i, j).RED() / 255.0);
            cc += col;
        }
        cc.clamp(0.0, 255.0);
        synthi.update(ii, jj, (hvColRGB<unsigned char>)cc);
    }
    sprintf(buff, "%s_synth.uhd.ppm", fname, ind);
    fd = fopen(buff, "wb");
    if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
    synthi.savePPM(fd, 1);
    fclose(fd);
    sprintf(buff, "%s_synth_wturb.uhd.ppm", fname, ind);
    fd = fopen(buff, "wb");
    if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
    synthib.savePPM(fd, 1);
    fclose(fd);

    sprintf(buff, "%s_wang_wturb.uhd.ppm", fname, ind);
    fd = fopen(buff, "wb");
    if (fd == 0) { perror("cannot open file:"); hvFatal(""); }
    synthiw.savePPM(fd, 1);
    fclose(fd);
    */

    return 0;
}


