//
// Created by grenier on 14/05/24.
//

#ifndef ASTEX_MESURE_PRESENCE_H
#define ASTEX_MESURE_PRESENCE_H

using namespace ASTex;

//-----------------------------------------------------------------------------
struct color_info{
    ImageRGBu8::PixelType couleur_;
    int compteur_;
    double compteurD_;

    color_info(ImageRGBu8::PixelType couleur){
        couleur_ = couleur;
        compteur_ = 1;
    }

    color_info(ImageRGBu8::PixelType couleur, int i){
        couleur_ = couleur;
        compteur_ = i;
    }

    color_info(ImageRGBu8::PixelType couleur, double i){
        couleur_ = couleur;
        compteurD_ = i;
    }

    color_info(ImageRGBu8::PixelType couleur, double i, int j){
        couleur_ = couleur;
        compteurD_ = i;
        compteur_ = j;
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

    bool operator==(const color_info& couleur2)
    {
        return couleur_ == couleur2.couleur_;
    }

    bool operator<(const color_info& couleur2)
    {
        return couleur_ < couleur2.couleur_;
    }
};


bool is_in(std::vector<color_info>& vecteur, ImageRGBu8::PixelType& couleur, int& id)
{
    id = -1;
//    for(ImageRGB8::PixelType col : vecteur)
    for(int i=0; i<vecteur.size(); i++)
    {
        if(vecteur.at(i) == couleur)
        {
            id = i;
            return true;
        }
    }
    return false;
}


bool is_in(std::vector<color_info>& vecteur, color_info& couleur, int& id)
{
    id = -1;
//    for(ImageRGB8::PixelType col : vecteur)
    for(int i=0; i<vecteur.size(); i++)
    {
        if(vecteur.at(i) == couleur.couleur_)
        {
            id = i;
            return true;
        }
    }
    return false;
}

void update_Tcontent(std::vector<color_info>& couleurs, color_info& col){
    int id = -1;
    bool presence = is_in(couleurs, col, id);

    if(not presence){
        couleurs.push_back(col);
    }
    else{
        couleurs.at(id).incr();
    }
}

void update_Tcontent(std::vector<color_info>& couleurs, color_info& col, int qtt){
    int id = -1;
    bool presence = is_in(couleurs, col, id);

    if(not presence){
        col.count(qtt);
        couleurs.push_back(col);
    }
    else{
        couleurs.at(id).incr(qtt);
    }
}

void update_Tcontent(std::vector<color_info>& couleurs, color_info& col, double qtt){
    int id = -1;
    bool presence = is_in(couleurs, col, id);

    if(not presence){
        col.count(qtt);
        couleurs.push_back(col);
    }
    else{
        couleurs.at(id).incr(qtt);
    }
}

#endif //ASTEX_MESURE_PRESENCE_H
