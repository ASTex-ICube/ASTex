
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



class pca {
public :
pca (const ImageRGBd& input){
    m_im.initItk(input.width(),input.height());
    m_im.itk() = input.itk();
}

void computeMeanColor(const mask_image& m){
    m_meancolor[0] = 0;
    m_meancolor[1] = 0;
    m_meancolor[2] = 0;
    int nbpix = 0;
    for (int x = 0; x < m_im.width(); ++x) {
        for (int y = 0; y < m_im.height(); ++y) {
            if (m.test(x,y))
            {
                ++nbpix;
                m_meancolor[0] += m_im.pixelAbsolute(x,y)[0];
                m_meancolor[1] += m_im.pixelAbsolute(x,y)[1];
                m_meancolor[2] += m_im.pixelAbsolute(x,y)[2];
            }
        }
    }
    m_meancolor[0] /= double(nbpix);
    m_meancolor[1] /= double(nbpix);
    m_meancolor[2] /= double(nbpix);
}

void computeCovariance(const mask_image& m){
    for (int z = 0; z < 3; ++z) {
        for (int t = 0; t < 3; ++t) {
            m_covar[z][t] = 0;
        }

    }

    int nbpix = 0;
    for (int x = 0; x < m_im.width(); ++x) {
        for (int y = 0; y < m_im.height(); ++y) {
            if (m.test(x,y))
            {
                ++nbpix;
                for (int z = 0; z < 3; ++z) {
                    for (int t = 0; t < 3; ++t) {
                        m_covar[z][t] += ( m_im.pixelAbsolute(x,y)[z] - m_meancolor[z] ) * ( m_im.pixelAbsolute(x,y)[t] - m_meancolor[t] );
                    }
                }
            }
        }
    }
    for (int z = 0; z < 3; ++z) {
        for (int t = 0; t < 3; ++t) {
            m_covar[z][t] /= double(nbpix);
        }
    }
}

void mat_mult(double vin[3], double vout[3]){
    for (int z = 0; z < 3; ++z) {
        vout[z] = m_covar[z][0]*vin[0] + m_covar[z][1]*vin[1] + m_covar[z][2]*vin[2];
    }
}

double dot(double v[3], double w[3]){
   return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

double normalize(double v[3]){
   double norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
   v[0] /= norm;
   v[1] /= norm;
   v[2] /= norm;
}

void cross(double v[3], double w[3], double res[3]){
    res[0] = v[1]*w[2] - v[2]*w[1];
    res[1] = v[2]*w[0] - v[0]*w[2];
    res[2] = v[0]*w[1] - v[1]*w[0];
}

void computeEigenVectors()
{
    /*** computing first EigenVector v1 = (a1 = theta, a2= phi in polar) ***/
    int i;
    double v [3];
    double w [3];
    double a1=0.0, a2=0.0;
    for(i=2;i<20;i++)
    {
        m_v1[0] = cos(a1); m_v1[1] = sin(a1)*cos(a2); m_v1[2] = sin(a1)*sin(a2); // cartesian coordinates
        mat_mult(m_v1,v); // matrix-vector mult by covar matrix
        w[0]=-sin(a1); w[1]=cos(a1)*cos(a2); w[2]=cos(a1)*sin(a2);
        double d1 = dot(v,w); // dot product
        w[0]=0.0; w[1]=-sin(a1)*sin(a2); w[2]=sin(a1)*cos(a2);
        double d2= dot(v,w); // dot product
        a1+=M_PI*(d1<0.0 ? (-1.0) : (1.0))/pow(2.0,(double)i);                /* adjust ai */
        a2+=M_PI*(d2<0.0 ? (-1.0) : (1.0))/pow(2.0,(double)i);
    }
    m_v1[0]=cos(a1); m_v1[1]=sin(a1)*cos(a2); m_v1[2]=sin(a1)*sin(a2); // first eigen vector

    /*** computing 2nd & 3rd EigenVectors (a1 = theta in polar in plane (v2,v3) ) ***/
    m_v2[0]=-m_v1[1]; m_v2[1]=m_v1[0]; m_v2[2]=0.0;
    normalize(m_v2);
    cross(m_v1,m_v2,m_v3);
    mat_mult(m_v2,v);
    double k1=-dot(v,m_v2);
    double k2=2.0*dot(v,m_v3);
    mat_mult(m_v3,v);
    k1 += dot(v,m_v3); /*   k1= v3.C.v3-v2.C.v2 ; k2= 2 * v3.C.v2   */
    a1=0.0; // ons single angle around v1
    for(i=2;i<20;i++)
        {
        double d1=sin(2.0*a1)*k1+cos(2.0*a1)*k2;
        a1 += M_PI*(d1<0.0 ? (-1.0) : (1.0))/pow(2.0,(double)i);
    }
    v[0]=m_v2[0]*cos(a1);v[1]=m_v2[1]*cos(a1);v[2]=m_v2[2]*cos(a1);
    w[0]=m_v3[0]*sin(a1);w[1]=m_v3[1]*sin(a1);w[2]=m_v3[2]*sin(a1);
    m_v2[0]=v[0]+w[0];m_v2[1]=v[1]+w[1];m_v2[2]=v[2]+w[2];
    normalize(m_v2);
    cross(m_v1,m_v2,m_v3);
//    hvMat3<T> res(v1,v2,v3);
//    res.transpose(); // TO DO : check with JMD if vectors are rows or columns
}

void computePCA (const mask_image& m){
    computeMeanColor(m);
    computeCovariance(m);
    computeEigenVectors();
}

void exportToTXT(std::string savefile){
    std::ofstream fichier(savefile, std::ios::out | std::ios::trunc);  // ouverture en écriture avec effacement du fichier ouvert
    fichier << "mean color"<< std::endl;
    fichier << m_meancolor[0] << "\t"<< m_meancolor[1] << "\t"<< m_meancolor[2] << std::endl;
    fichier << "eigen vectors"<< std::endl;
    fichier << m_v1[0] << "\t"<< m_v1[1] << "\t"<< m_v1[2] << std::endl;
    fichier << m_v2[0] << "\t"<< m_v2[1] << "\t"<< m_v2[2] << std::endl;
    fichier << m_v3[0] << "\t"<< m_v3[1] << "\t"<< m_v3[2] << std::endl;
    fichier.close();
}

void project(ImageGrayd& res1, ImageGrayd& res2, ImageGrayd& res3){
    assert(m_im.width() == res1.width() && m_im.height() == res1.height());
    assert(m_im.width() == res2.width() && m_im.height() == res2.height());
    assert(m_im.width() == res3.width() && m_im.height() == res3.height());

    for (int x = 0; x < m_im.width(); ++x) {
        for (int y = 0; y < m_im.height(); ++y) {
            res1.pixelAbsolute(x,y) = (m_im.pixelAbsolute(x,y)[0]-m_meancolor[0]) * m_v1[0] + (m_im.pixelAbsolute(x,y)[1]-m_meancolor[1]) * m_v1[1] + (m_im.pixelAbsolute(x,y)[2]-m_meancolor[2]) * m_v1[2];
            res2.pixelAbsolute(x,y) = (m_im.pixelAbsolute(x,y)[0]-m_meancolor[0]) * m_v2[0] + (m_im.pixelAbsolute(x,y)[1]-m_meancolor[1]) * m_v2[1] + (m_im.pixelAbsolute(x,y)[2]-m_meancolor[2]) * m_v2[2];
            res3.pixelAbsolute(x,y) = (m_im.pixelAbsolute(x,y)[0]-m_meancolor[0]) * m_v3[0] + (m_im.pixelAbsolute(x,y)[1]-m_meancolor[1]) * m_v3[1] + (m_im.pixelAbsolute(x,y)[2]-m_meancolor[2]) * m_v3[2];
        }
    }
}

// image
ImageRGBd m_im;
    // eigenvectors
double m_v1 [3];
double m_v2 [3];
double m_v3 [3];
  // mean values
double m_meancolor [3];
// covariance matrix
double m_covar [3][3];
};

/**
 * @brief RPnoise random phase noise
 * @param inputfile a grayscale input example
 * @param outputfile is a patchwork image containing input + spectrum + output (pure random phase)
 * @return error code
 */
int ACP(std::string name_file, std::string inputfile, int size_fft, float sig_freq,  std::vector<ImageGrayd> masks)
{
    std::cout << "TRAITEMENT DE  " <<inputfile<< std::endl;


    // LOAD INPUT

    //********** COULEUR ***********/

    ImageRGBd input_color;

    input_color.load(inputfile);

    typedef ImageRGBd::ItkImg IMG_DBL;
    typedef ImageRGBf::ItkImg IMG_FLT;
    typedef ImageRGBu8::ItkImg IMG_U8;

    // first transform [0,255] double -> [0,1] double
    ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::Pointer filter0 =
            ColorSpace::FilterRGB255To01<IMG_DBL,IMG_DBL>::New();
    filter0->SetInput(input_color.itk());
    filter0->Update();

    input_color.itk() = filter0->GetOutput();



    //OUTPUT PCA
    ImageGrayd coord1;
    coord1.initItk(input_color.width(),input_color.height());
    ImageGrayd coord2;
    coord2.initItk(input_color.width(),input_color.height());
    ImageGrayd coord3;
    coord3.initItk(input_color.width(),input_color.height());
    /************* Fichier de save ******************/


    //Creation file tree
    std::string path = "/home/guingo/Seafile/Geoffrey/ACP/";
    system(std::string("mkdir -p "+path).c_str());
    path = "/home/guingo/Seafile/Geoffrey/ACP/"+name_file+"/"+std::to_string(size_fft)+"/"+std::to_string(sig_freq)+"/";
    system(std::string("mkdir -p "+path).c_str());



    /********************* PCA ************************/
    pca pca(input_color); // Input : sensé être notre noise RGB

    std::vector<ImageGrayd> masks_normalized;
    masks_normalized.resize(masks.size());

    normalize_distance_maps(masks,masks_normalized);

            //Pour chaque masque
    for (int k = 0; k < masks.size(); ++k) {

        //Creation de notre masque et export du binaire
        mask_largest_values m(masks_normalized[k],0.3);

        ImageGrayu8 binary;
        binary.initItk(input_color.width(),input_color.height());
        m.export_binary_image(binary);
        binary.save(path+name_file+"_nclusters_"+std::to_string(masks.size())+"_binary_mask_"+std::to_string(k)+".png");

        //Calcul de la base
        pca.computePCA(m);
        //Export de la base
        std::string savefile = path+name_file+"_nclusters_"+
        std::to_string(masks.size())+ "_fftsize_"+std::to_string(size_fft)+"_ACP_mask_"+std::to_string(k)+".txt";

        pca.exportToTXT(savefile);

        pca.project(coord1,coord2,coord3);

        std::vector<ImageGrayd> spectrum;
        spectrum.resize(3);

        // Charger les spectres par auto-covariance

        spectrum[0].initItk(64,64);
        spectrum[1].initItk(64,64);
        spectrum[2].initItk(64,64);

        spectrum_by_autocorrelation_small_size(coord1,spectrum[0],m,0);
        spectrum_by_autocorrelation_small_size(coord2,spectrum[1],m,0);
        spectrum_by_autocorrelation_small_size(coord3,spectrum[2],m,0);

        save_spectrum(spectrum[0],path+name_file+"_spectrum_nclusters_"+
                std::to_string(masks.size())+"_mask_"+std::to_string(k)+"_coord0.png");
        save_spectrum(spectrum[1],path+name_file+"_spectrum_nclusters_"+
                std::to_string(masks.size())+"_mask_"+std::to_string(k)+"_coord1.png");
        save_spectrum(spectrum[2],path+name_file+"_spectrum_nclusters_"+
                std::to_string(masks.size())+"_mask_"+std::to_string(k)+"_coord2.png");

    }




     /*********************************************/

    return EXIT_SUCCESS;
}

int main( int argc, char ** argv )
{
    if( argc < 7 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " <Source name> <Source path>  <Size_FFT> <SIGfilter> <MAsk1> <Mask2> [<MaskN>]*" << std::endl;
    return EXIT_FAILURE;
    }

    const int nb_maps = (argc - 5) ;

    std::string filename = argv[1];
    std::string path_source = argv[2];

    int size_fft = atoi(argv[3]);
    float sig_freq = atof(argv[4]);

    std::vector<ImageGrayd> masks;
    masks.resize(nb_maps);

    //ATTENTION : Les masques doivent être des masque de coordonnées positives pour fit le load
    for (int i=0; i<nb_maps; ++i)
       load(masks[i],argv[i+5],0,2);

    std::cout << " ok go " << filename << " " << path_source << std ::endl;

    return ACP(filename,path_source,size_fft,sig_freq,masks);
}
