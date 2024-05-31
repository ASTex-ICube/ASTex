//
// Created by grenier on 23/05/24.
//

#ifndef ASTEX_MESURE_TRIPLET_H
#define ASTEX_MESURE_TRIPLET_H

using namespace ASTex;

struct color_triple{
    ImageRGBu8::PixelType couleur1_;
    ImageRGBu8::PixelType couleur2_;
    ImageRGBu8::PixelType couleur3_;
    int compteur_;
    double compteurD_;

    color_triple();

    color_triple(ImageRGBu8::PixelType couleur1, ImageRGBu8::PixelType couleur2, ImageRGBu8::PixelType couleur3){
        couleur1_ = couleur1;
        couleur2_ = couleur2;
        couleur3_ = couleur3;
        compteur_ = 1;
    }

    color_triple(ImageRGBu8::PixelType couleur1, ImageRGBu8::PixelType couleur2, ImageRGBu8::PixelType couleur3, int i){
        couleur1_ = couleur1;
        couleur2_ = couleur2;
        couleur3_ = couleur3;
        compteur_ = i;
    }

    color_triple(ImageRGBu8::PixelType couleur1, ImageRGBu8::PixelType couleur2, ImageRGBu8::PixelType couleur3, double i){
        couleur1_ = couleur1;
        couleur2_ = couleur2;
        couleur3_ = couleur3;
        compteurD_ = i;
    }

    void incr()
    {
        compteur_ += 1;
    }

    void incr(int i)
    {
        compteur_ += i;
    }

    void incr(double i)
    {
        compteurD_ += i;
    }

    void count(int qtt){
        compteur_ = qtt;
    }

    void count(double qtt){
        compteurD_ = qtt;
    }

    bool operator==(const color_triple& couleur2)
    {
        return (couleur1_ == couleur2.couleur1_) and (couleur2_ == couleur2.couleur2_) and (couleur3_ == couleur2.couleur3_);
    }

    bool operator<(const color_triple& couleur2)
    {
        if(couleur1_ == couleur2.couleur1_ and couleur2_ == couleur2.couleur2_){
            return (couleur3_ < couleur2.couleur3_);
        }
        else if(couleur1_ == couleur2.couleur1_){
            return (couleur2_ < couleur2.couleur2_);
        }
        else{
            return (couleur1_ < couleur2.couleur1_);
        }
    }
};

bool is_in(std::vector<color_triple>& vecteur, ImageRGBu8::PixelType& couleur1, ImageRGBu8::PixelType& couleur2, ImageRGBu8::PixelType& couleur3, int& id)
{
    id = -1;
//    for(ImageRGB8::PixelType col : vecteur)
    for(int i=0; i<vecteur.size(); i++)
    {
        if(vecteur.at(i) == color_triple{couleur1, couleur2, couleur3})
        {
            id = i;
            return true;
        }
    }
    return false;
}

bool is_in(std::vector<color_triple>& vecteur, color_triple& couleurs, int& id)
{
    id = -1;
//    for(ImageRGB8::PixelType col : vecteur)
    for(int i=0; i<vecteur.size(); i++)
    {
        if(vecteur.at(i) == couleurs)
        {
            id = i;
            return true;
        }
    }
    return false;
}

bool get_triplet(color_triple& triplet, ImageRGBu8::PixelType& couleur1, ImageRGBu8::PixelType& couleur2, ImageRGBu8::PixelType& couleur3, ImageRGBu8::PixelType& couleur4){
    std::vector couleurs{couleur1, couleur2, couleur3, couleur4};
    std::sort(couleurs.begin(), couleurs.end());

    auto last = std::unique(couleurs.begin(), couleurs.end());
    couleurs.erase(last, couleurs.end());

    if(couleurs.size()!=3){return false;}
    else{
        triplet = color_triple{couleurs.at(0), couleurs.at(1), couleurs.at(2)};
        return true;
    }
}

void update_Tcontent(std::vector<color_triple>& couleurs, color_triple& triplet){
    int id = -1;
    bool presence = is_in(couleurs, triplet, id);

    if(not presence){
        couleurs.push_back(triplet);
    }
    else{
        couleurs.at(id).incr();
    }
}

void update_Tcontent(std::vector<color_triple>& couleurs, color_triple& triplet, int qtt){
    int id = -1;
    bool presence = is_in(couleurs, triplet, id);

    if(not presence){
        triplet.count(qtt);
        couleurs.push_back(triplet);
    }
    else{
        couleurs.at(id).incr(qtt);
    }
}


void update_Tcontent(std::vector<color_triple>& couleurs, color_triple& triplet, double qtt){
    int id = -1;
    bool presence = is_in(couleurs, triplet, id);

    if(not presence){
        triplet.count(qtt);
        couleurs.push_back(triplet);
    }
    else{
        couleurs.at(id).incr(qtt);
    }
}


#endif //ASTEX_MESURE_TRIPLET_H
