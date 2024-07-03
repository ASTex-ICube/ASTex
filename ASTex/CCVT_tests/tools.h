//
// Created by grenier on 26/06/23.
//

#ifndef ASTEX_TOOLS_H
#define ASTEX_TOOLS_H

#include <ASTex/easy_io.h>
#include <ASTex/image_rgb.h>


using namespace ASTex;

void save_pgm(const std::string& filename, ImageGrayd& data)
{
    std::ofstream output(filename.c_str());
    output << "P2" << std::endl;
    output << "# " << filename << std::endl;
    output << data.width() << " " << data.height() << std::endl;
    output << "255" << std::endl;
    data.for_all_pixels([&output](typename ImageGrayd::PixelType& P, int x, int y)
                        {
                            int value = static_cast<int>(255*P);
                            output << value << std::endl;
                        });
    output.close();
}






// ---------------------------------------------------------------------------
double gauss(double moy1, double moy2, double var1, double var2, double x, double y, double scale)
{

    double mu1 = scale*moy1;
    double mu2 = scale*moy2;
    double sig1 = scale*scale*var1;
    double sig2 = scale*scale*var2;

    // gaussienne en x
    double X = x-mu1;
    double G1 = std::exp(-(X*X)/(2.*sig1))/(std::sqrt(2.*M_PI*sig1));

    // gaussienne en y
    double Y = y-mu2;
    double G2 = std::exp(-(Y*Y)/(2.*sig2))/(std::sqrt(2.*M_PI*sig2));

    return G1*G2;
}



double gaussAC_2D(double moy, double var, double AC, double x, double y)
{
    double ecart_type = std::sqrt(var);
    /* matrice covariance :
     * var AC
     * AC  var
     */
    double det = (var*var)-(AC*AC); // déterminant de la matrice de covariance

    /* matrice précision (inverse matrice covariance) :
     * var/det -AC/det
     * -AC/det  var/det
     */
    double var_inv = var/det;
    double AC_inv = -AC/det;


    double X = x-moy;
    double Y = y-moy;

    double produit = X*(X*var_inv + Y*AC_inv) + Y*(X*AC_inv + Y*var_inv);
    double norm = 2.*M_PI*std::sqrt(det);

    double G = std::exp(-0.5*produit)/norm;

    return G;
}





struct Graine{
    double x_;
    double y_;
    double weight_;

    Graine(){
        x_=0.;
        y_=0.;
        weight_ = 1.;
    }

    Graine(double x, double y, double w){
        x_ = x;
        y_ = y;
        weight_ = w;
    }
};


double grain_dist(double couleur1, double couleur2, std::vector<double> H_color, std::vector<Graine> H_seeds){
    int id1 = -1;
    int id2 = -1;
    for(int i=0; i<H_color.size(); i++){
        if(H_color.at(i)==couleur1){id1=i;}
        if(H_color.at(i)==couleur2){id2=i;}
    }

    Graine g1 = H_seeds.at(id1);
    Graine g2 = H_seeds.at(id2);

    return sqrt((g1.x_-g2.x_)*(g1.x_-g2.x_) + (g1.y_-g2.y_)*(g1.y_-g2.y_));
}






#endif //ASTEX_TOOLS_H
