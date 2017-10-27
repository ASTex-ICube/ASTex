/*******************************************************************************
* ASTex:                                                                       *
* Copyright (C) IGG Group, ICube, University of Strasbourg, France             *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: https://astex-icube.github.io                                      *
* Contact information: astex@icube.unistra.fr                                  *
*                                                                              *
*******************************************************************************/



#include "io.h"
#include "fourier.h"
#include "itkJPEGImageIO.h"
#include "mode_seeking.h"
#include "output_geoffrey.h"
#include "image_treatment.h"
#include "utilities.h"
#include "tests.h"
#include "colorspace_filters.h"


using namespace ASTex;


int spectrum_reset (ImageGrayd& s)
{
    for (int fx=0; fx < s.width(); ++fx)
    {
        for (int fy=0; fy < s.height(); ++fy)
        {
            s.pixelAbsolute(fx,fy) = 0.0;
        }
    }
    return 0;
}

int spectrum_accumulate (ImageGrayd& accu, const ImageGrayd& s)
{
    for (int fx=0; fx < s.width(); ++fx)
    {
        for (int fy=0; fy < s.height(); ++fy)
        {
            accu.pixelAbsolute(fx,fy) += s.pixelAbsolute(fx,fy);
        }
    }
    return 0;
}

int spectrum_rescale (ImageGrayd& s, double ratio)
{
    for (int fx=0; fx < s.width(); ++fx)
    {
        for (int fy=0; fy < s.height(); ++fy)
        {
            s.pixelAbsolute(fx,fy) *= ratio;
        }
    }
    return 0;
}

int distance_invert_expo (const ImageGrayd& s, ImageGrayd& d, double dev)
{
    for (int fx=0; fx < s.width(); ++fx)
    {
        for (int fy=0; fy < s.height(); ++fy)
        {
            d.pixelAbsolute(fx,fy) = std::exp(-s.pixelAbsolute(fx,fy)/dev);
        }
    }
    return 0;
}

int Image_log (const ImageGrayd& s, ImageGrayd& d, double ratio)
{
    for (int fx=0; fx < s.width(); ++fx)
    {
        for (int fy=0; fy < s.height(); ++fy)
        {
            d.pixelAbsolute(fx,fy) = std::sqrt((s.pixelAbsolute(fx,fy)*ratio )) ;
        }
    }
    return 0;
}

int Image_expo (const ImageGrayd& s, ImageGrayd& d, double dev)
{
    for (int fx=0; fx < s.width(); ++fx)
    {
        for (int fy=0; fy < s.height(); ++fy)
        {
            d.pixelAbsolute(fx,fy) = 1- std::exp(-s.pixelAbsolute(fx,fy)/dev);
        }
    }
    return 0;
}

double get_average(const ImageGrayd s)
{
    double sum=0;
    for (int fx=0; fx < s.width(); ++fx)
    {
        for (int fy=0; fy < s.height(); ++fy){
            sum += s.pixelAbsolute(fx,fy);
        }
    }

    return sum/(s.width()*s.height());

}

double get_max(const ImageGrayd s)
{
    double max=0;
    for (int fx=0; fx < s.width(); ++fx)
    {
        for (int fy=0; fy < s.height(); ++fy){
            if( max <s.pixelAbsolute(fx,fy) )
                max = s.pixelAbsolute(fx,fy);
        }
    }

    return max;
}



/**
 * @brief kmean
 * @param input a grayscale input example
 * @param
 * @return error code
 */
int kmean(ImageGrayd input, std::vector<ImageGrayd>& output_masks, std::vector<ImageGrayd>& masks_spectrum_save, std::vector<ImageGrayd>& masks_spectrum_visu, const int nb_clusters,const int welchsize = 15,const int welchstep = 4,const int fft_size=8, std::string name_file="test")
{
    srand(time(NULL));
    const int im_w = input.width();
    const int im_h = input.height();

    for (int c=0; c<nb_clusters; ++c)
    {
        output_masks[c].initItk(im_w,im_h,true);
        masks_spectrum_visu[c].initItk(fft_size,fft_size,true);
        masks_spectrum_save[c].initItk(fft_size,fft_size,true);

    }



    std::cout << "compute local spectrum ...\r" << std::flush;
    local_spectrum lsp (input,fft_size,'p');
    lsp.welch_post_loading(welchsize, welchstep);
    std::cout << "compute local spectrum - done" << std::endl;

    ImageGray16 clusters; // image of cluster indices
    clusters.initItk(im_w, im_h, true);

    std::vector< ImageGrayd > seeds; // seeds = spectra
    seeds.resize(nb_clusters);

    std::vector<int> clusters_size; // nb elements in each cluster
    clusters_size.resize(nb_clusters);

    std::vector< ImageGrayd > seeds_distance_maps; 	// allocate distance maps
    seeds_distance_maps.resize(nb_clusters);

    std::vector< std::vector<double> > distance_inter_seed; 	//
    distance_inter_seed.resize( nb_clusters , std::vector<double>( nb_clusters ) );

    std::vector< double > distance_intra_classe; 	//
    distance_intra_classe.resize(nb_clusters);




    std::vector< double > average_dist_intra_class; 	//
    average_dist_intra_class.resize(nb_clusters);

    ImageGrayd sp_tmp; // tempo spectrum
    sp_tmp.initItk(fft_size,fft_size,true);

    for (int c = 0; c < nb_clusters; ++c)
    {
        seeds[c].initItk(fft_size,fft_size,true);
        clusters_size[c] = -1;
        seeds_distance_maps[c].initItk(im_w,im_h,true);
    }

    // random seed init
    for (int c = 0; c < nb_clusters; ++c)
    {
        int x = rand()% im_w;
        int y = rand()% im_h;
        lsp.get_spectrum(seeds[c],x,y);
    }

    // iterate

    int iter = 0;
    const int max_iter = 20;
    int nb_changes = 1;

    while (iter < max_iter && nb_changes > 0)
    {
        ++iter;
        nb_changes = 0;
        std::cout << "iteration " << iter << " ...\r" << std::flush;

        // compute distance maps
        // std::cout << "iteration " << iter << " ..." << std::endl;
        // std::cout << "compute distance maps" << std::endl;
        for (int c = 0; c < nb_clusters; ++c)
        {
//            lsp.distance_map_to_spectum_linear_weights(seeds_distance_maps[c],seeds[c]);
            lsp.distance_map_to_spectum_uniform_weights(seeds_distance_maps[c],seeds[c]);

        }


        // compute clusters
        // std::cout << "compute clusters" << std::endl;
        for (int x=0; x<im_w; ++x)
        {
            for (int y=0; y<im_h; ++y)
            {
                int c_min = 0;
                double dist_min = seeds_distance_maps[0].pixelAbsolute(x,y);
                for (int c = 1; c < nb_clusters; ++c)
                {
                    double dist = seeds_distance_maps[c].pixelAbsolute(x,y);
                    if (dist < dist_min)
                    {
                        c_min = c;
                        dist_min = dist;
                    }
                }
                if (clusters.pixelAbsolute(x,y) != c_min)
                {
                    nb_changes++;
                }
                clusters.pixelAbsolute(x,y) = c_min;
            }
        }

        // update seeds
        // std::cout << "reset seeds" << std::endl;
        for (int c = 0; c < nb_clusters; ++c)
        {
            // std::cout << c << std::endl;
            spectrum_reset(seeds[c]);
            clusters_size[c] = 0;
        }

        // std::cout << "accu seeds" << std::endl;
        for (int x=0; x<im_w; ++x)
        {
            for (int y=0; y<im_h; ++y)
            {
                // std::cout << "x= "<< x <<" / y= " << y << std::endl;

                const int c = clusters.pixelAbsolute(x,y);
                lsp.get_spectrum(sp_tmp,x,y);
                spectrum_accumulate(seeds[c],sp_tmp);
                clusters_size[c] ++;
            }
        }

        // std::cout << "rescale seeds" << std::endl;
        for (int c = 0; c < nb_clusters; ++c)
        {
            if (clusters_size[c] != 0)
            {
                spectrum_rescale(seeds[c], 1.0 / double(clusters_size[c]));
            }
            else
            {
                std::cout << "cluster " << c << " is empty" << std::endl;
            }
        }

        std::cout << "iteration " << iter << " - done - "<< nb_changes << " pixels changed" << std::endl;

    }
     // KMEAN TERMINE
    for (int c=0; c<nb_clusters; ++c)
    {
        output_masks[c].initItk(im_w,im_h,true);
        masks_spectrum_visu[c].initItk(fft_size,fft_size,true);
        masks_spectrum_save[c].initItk(fft_size,fft_size,true);

    }

    std::string savefile = "/home/guingo/Images/Test_methods/bench_article/Kmean/"+name_file+"/"+name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size)+"_COORD_DIST_uniformkmean.txt";
    std::ofstream fichier(savefile, std::ios::out | std::ios::trunc);  // ouverture en écriture avec effacement du fichier ouvert

    std::string savefile_positive = "/home/guingo/Images/Test_methods/bench_article/Kmean/"+name_file+"/"+name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size)+"_COORD_positive_DIST_uniformkmean.txt";
    std::ofstream fichier_positive(savefile_positive, std::ios::out | std::ios::trunc);  // ouverture en écriture avec effacement du fichier ouvert


    std::string savefile_uniforme = "/home/guingo/Images/Test_methods/bench_article/Kmean/"+name_file+"/"+name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size)+"_COORD_uniforme_DIST_uniformkmean.txt";
    std::ofstream fichier_uniforme(savefile_uniforme, std::ios::out | std::ios::trunc);  // ouverture en écriture avec effacement du fichier ouvert


    std::string savefile_positive_uniform = "/home/guingo/Images/Test_methods/bench_article/Kmean/"+name_file+"/"+name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size)+"_COORD_positive_uniforme_DIST_uniformkmean.txt";
    std::ofstream fichier_positive_uniforme(savefile_positive_uniform, std::ios::out | std::ios::trunc);  // ouverture en écriture avec effacement du fichier ouvert

    std::vector<ImageGrayd> masks_coord;
    masks_coord.resize(nb_clusters);

    std::vector<ImageGrayd> masks_coord_uniform;
    masks_coord_uniform.resize(nb_clusters);

    std::vector<ImageGrayd> masks_coord_positive;
    masks_coord_positive.resize(nb_clusters);

    std::vector<ImageGrayd> masks_coord_positive_uniform;
    masks_coord_positive_uniform.resize(nb_clusters);

    std::vector<ImageGrayd> masks_coord_load;
    masks_coord_load.resize(nb_clusters);

    std::vector<ImageGrayd> masks_coord_substract;
    masks_coord_substract.resize(nb_clusters);

    if(nb_clusters == 2){
        spectrum_projector_2D proj2D;

        for (int c = 0; c < nb_clusters; ++c){
            proj2D.add(seeds[c]);
             masks_coord[c].initItk(im_w,im_h,true);
             masks_coord_substract[c].initItk(im_w,im_h,true);
             masks_coord_positive[c].initItk(im_w,im_h,true);
             masks_coord_positive_uniform[c].initItk(im_w,im_h,true);
             masks_coord_uniform[c].initItk(im_w,im_h,true);
        }

        proj2D.precompute_system_uniform_weights();

        for (int x=0; x<im_w; ++x){
            for (int y=0; y<im_h; ++y){

                std::vector<double>coord;

                ImageGrayd sp_tmp; // tempo spectrum
                sp_tmp.initItk(fft_size,fft_size,true);

                lsp.get_spectrum(sp_tmp,x,y);

                proj2D.project_uniform_weights(sp_tmp,coord);
                for (int c = 0; c < nb_clusters; ++c){
                    masks_coord_uniform[c].pixelAbsolute(x,y)=coord[c];
                 }
                fichier_uniforme << x << "\t" << y <<"\t"<<clusters.pixelAbsolute(x,y)<<"\t"<< coord[0]<<"\t"<< coord[1]<<std::endl;

                proj2D.reproject_coordinates_to_positive(coord);

                for (int c = 0; c < nb_clusters; ++c){
                    masks_coord_positive_uniform[c].pixelAbsolute(x,y)=coord[c];
                }
                fichier_positive_uniforme << x << "\t" << y <<"\t"<<clusters.pixelAbsolute(x,y)<<"\t"<< coord[0]<<"\t"<< coord[1]<<std::endl;
            }
        }
        proj2D.precompute_system();

        for (int x=0; x<im_w; ++x){
            for (int y=0; y<im_h; ++y){

                std::vector<double>coord;

                ImageGrayd sp_tmp; // tempo spectrum
                sp_tmp.initItk(fft_size,fft_size,true);

                lsp.get_spectrum(sp_tmp,x,y);

                proj2D.project(sp_tmp,coord);
                for (int c = 0; c < nb_clusters; ++c){
                    masks_coord[c].pixelAbsolute(x,y)=coord[c];
                 }
                fichier << x << "\t" << y <<"\t"<<clusters.pixelAbsolute(x,y)<<"\t"<< coord[0]<<"\t"<< coord[1]<<std::endl;

                proj2D.reproject_coordinates_to_positive(coord);

                for (int c = 0; c < nb_clusters; ++c){
                    masks_coord_positive[c].pixelAbsolute(x,y)=coord[c];
                }
                fichier_positive << x << "\t" << y <<"\t"<<clusters.pixelAbsolute(x,y)<<"\t"<< coord[0]<<"\t"<< coord[1]<<std::endl;
            }
        }
    }else if (nb_clusters ==3) {
    spectrum_projector_3D proj3D;

    for (int c = 0; c < nb_clusters; ++c){
        proj3D.add(seeds[c]);
         masks_coord[c].initItk(im_w,im_h,true);
         masks_coord_substract[c].initItk(im_w,im_h,true);
         masks_coord_positive[c].initItk(im_w,im_h,true);
         masks_coord_positive_uniform[c].initItk(im_w,im_h,true);
         masks_coord_uniform[c].initItk(im_w,im_h,true);
    }

    proj3D.precompute_system_uniform_weights();

    for (int x=0; x<im_w; ++x){
        for (int y=0; y<im_h; ++y){

            std::vector<double>coord;

            ImageGrayd sp_tmp; // tempo spectrum
            sp_tmp.initItk(fft_size,fft_size,true);

            lsp.get_spectrum(sp_tmp,x,y);

            proj3D.project_uniform_weights(sp_tmp,coord);
            for (int c = 0; c < nb_clusters; ++c){
                masks_coord_uniform[c].pixelAbsolute(x,y)=coord[c];
             }
            fichier_uniforme << x << "\t" << y <<"\t"<<clusters.pixelAbsolute(x,y)<<"\t"<< coord[0]<<"\t"<< coord[1]<<"\t"<< coord[2]<<std::endl;
            proj3D.reproject_coordinates_to_positive(coord);

            for (int c = 0; c < nb_clusters; ++c){
                masks_coord_positive_uniform[c].pixelAbsolute(x,y)=coord[c];
            }
            fichier_positive_uniforme << x << "\t" << y <<"\t"<<clusters.pixelAbsolute(x,y)<<"\t"<< coord[0]<<"\t"<< coord[1]<<"\t"<< coord[2]<<std::endl;

        }
    }

    proj3D.precompute_system();

    for (int x=0; x<im_w; ++x){
        for (int y=0; y<im_h; ++y){

            std::vector<double>coord;

            ImageGrayd sp_tmp; // tempo spectrum
            sp_tmp.initItk(fft_size,fft_size,true);

            lsp.get_spectrum(sp_tmp,x,y);
            proj3D.project(sp_tmp,coord);

            for (int c = 0; c < nb_clusters; ++c){
                masks_coord[c].pixelAbsolute(x,y)=coord[c];
            }
            fichier << x << "\t" << y <<"\t"<<clusters.pixelAbsolute(x,y)<<"\t"<< coord[0]<<"\t"<< coord[1]<<"\t"<< coord[2]<<std::endl;

            proj3D.reproject_coordinates_to_positive(coord);

            for (int c = 0; c < nb_clusters; ++c){
                masks_coord_positive[c].pixelAbsolute(x,y)=coord[c];
            }
            fichier_positive << x << "\t" << y <<"\t"<<clusters.pixelAbsolute(x,y)<<"\t"<< coord[0]<<"\t"<< coord[1]<<"\t"<< coord[2]<<std::endl;



        }
    }
    }


    std::string path_out="/home/guingo/Images/Test_methods/bench_article/Kmean/"+name_file+"/distance_coords/";
    system(std::string("mkdir -p "+path_out).c_str());

    if(nb_clusters ==2 ){
    for (int c=0; c<nb_clusters; ++c){
            save(masks_coord[c], path_out+name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size) + "_coord_"+std::to_string(c)+".png",-1,3);
            save(masks_coord_positive[c], path_out+name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size) + "_coord_positive_"+std::to_string(c)+".png",0,1);
            save(masks_coord_uniform[c], path_out+name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size) + "_coord_"+std::to_string(c)+".png",-1,3);
            save(masks_coord_positive_uniform[c], path_out+name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size) + "_coord_positive_"+std::to_string(c)+".png",0,1);

        }
    }else{
        for (int c=0; c<nb_clusters; ++c){
                save(masks_coord[c], path_out+name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size) + "_coord_"+std::to_string(c)+".png",-1,3);
                save(masks_coord_positive[c], path_out+name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size) + "_coord_positive_"+std::to_string(c)+".png",0,1);
                save(masks_coord_uniform[c], path_out+name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size) + "_coord_"+std::to_string(c)+".png",-1,3);
                save(masks_coord_positive_uniform[c], path_out+name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size) + "_coord_positive_"+std::to_string(c)+".png",0,1);
            }

    }

//        for (int c=0; c<nb_clusters; ++c){
//            load(masks_coord_load[c], path_out+name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size) + "_coord_positive_"+std::to_string(c)+".png",-1,3);
////            save_coordinates(masks_coord_load[c], path_out+name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size) + "_coord2_"+std::to_string(c)+".png");
//      }


//       for (int c=0; c<nb_clusters; ++c){

//           for (int i = 0; i < im_w; ++i)
//               for (int j = 0; j < im_h; ++j)
//                   masks_coord_substract[c].pixelAbsolute(i,j) = masks_coord_load[c].pixelAbsolute(i,j)- masks_coord[c].pixelAbsolute(i,j);

//       //ADD 127 TO substract
//           for (int i = 0; i < im_w; ++i)
//               for (int j = 0; j < im_h; ++j)
//                   masks_coord_substract[c].pixelAbsolute(i,j) +=clamp( masks_coord_substract[c].pixelAbsolute(i,j)+0.5,0,1);

//           print_image(masks_coord_substract[c],path_out,name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size) + "_coord_substract_"+std::to_string(c), 1.0);

//       }

//       for (int c=0; c<nb_clusters; ++c){
//           load_coordinates(masks_coord[c], path_out+name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size) + "_coord2_"+std::to_string(c)+".png");
//       }

//       for (int c=0; c<nb_clusters; ++c){

//           for (int i = 0; i < im_w; ++i)
//               for (int j = 0; j < im_h; ++j)
//                   masks_coord_substract[c].pixelAbsolute(i,j) = masks_coord_load[c].pixelAbsolute(i,j)- masks_coord[c].pixelAbsolute(i,j);

//       //ADD 127 TO substract
//           for (int i = 0; i < im_w; ++i)
//               for (int j = 0; j < im_h; ++j)
//                   masks_coord_substract[c].pixelAbsolute(i,j) +=clamp( masks_coord_substract[c].pixelAbsolute(i,j)+0.5,0,1);

//           print_image(masks_coord_substract[c],path_out,name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size) + "_coord2_substract_"+std::to_string(c), 1.0);

//       }


//        print_image(masks_coord_load[c],path_out,name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size) + "_coord_"+std::to_string(c), 1.0);

    fichier.close();
    fichier_positive.close();

    return 0;
}




int main( int argc, char ** argv )
{
    if( argc < 3 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " <Source> <nb clusters>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string filename_source = argv[1];
    const int nb_clusters = atoi (argv[2]);

    std::string delimiter = ".p";
    std::string name_file = filename_source.substr(0, filename_source.find(delimiter));

    std::string path_in_b ="/home/guingo/Images/Test_methods/bench_article/Kmean/"+name_file+"/";
    std::string path_out="/home/guingo/Images/Test_methods/bench_article/Kmean/"+name_file+"/distance_maps/";
    std::string path_norm="/home/guingo/Images/Test_methods/bench_article/Kmean/"+name_file+"/kmean_mask/";
    std::string path_spectrum_visu="/home/guingo/Images/Test_methods/bench_article/Kmean/"+name_file+"/spectrum_visu/";
    std::string path_spectrum_save="/home/guingo/Images/Test_methods/bench_article/Kmean/"+name_file+"/spectrum_save/";


    system(std::string("mkdir -p "+path_in_b).c_str());
    system(std::string("mkdir -p "+path_out).c_str());
    system(std::string("mkdir -p "+path_norm).c_str());
    system(std::string("mkdir -p "+path_spectrum_visu).c_str());
    system(std::string("mkdir -p "+path_spectrum_save).c_str());


    std::cout << "TRAITEMENT "<<name_file << std::endl;

    // LOAD INPUT
    ImageGrayd input;
    load2lightness(filename_source, input);

    int fft_size = 2;

    for (fft_size = 2; fft_size <= 32; fft_size*=2)
    {
        std::cout << name_file << "fft_size  "<< fft_size<<" NB CLUSTER "<<nb_clusters<< std::endl;

        // COMPUTE MASKS
        std::vector<ImageGrayd> masks;
        masks.resize(nb_clusters);
        std::vector<ImageGrayd> masks_norm;
        masks_norm.resize(nb_clusters);
        std::vector<ImageGrayd> masks_spectrum_save;
        masks_spectrum_save.resize(nb_clusters);
        std::vector<ImageGrayd> masks_spectrum_visu;
        masks_spectrum_visu.resize(nb_clusters);

        kmean(input,masks,masks_spectrum_save,masks_spectrum_visu,nb_clusters,(fft_size*2)-1,4,fft_size, name_file );

        // SAVE MASKS
        // std::string path_in_b ="/home/guingo/Images/Test_methods/bench_article/";
//        normalize_distance_maps(masks, masks_norm);

//        for (int c=0; c<nb_clusters; ++c)
//        {
//            print_image(masks[c],path_out,name_file+"_nclusters_"+ std::to_string(nb_clusters)+ "_fftsize_"+std::to_string(fft_size) + "_mask_"+std::to_string(c), 1.0);

//        }


    }


    return 0;
}
