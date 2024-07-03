//
// Created by grenier on 14/05/24.
//

#ifndef ASTEX_MESURE_VOISINAGE_H
#define ASTEX_MESURE_VOISINAGE_H

using namespace ASTex;

// ---------------------------------------------------------------------------
struct color_vois{
    ImageRGBu8::PixelType couleur1_;
    ImageRGBu8::PixelType couleur2_;
    int compteur_;
    double compteurD_;

    color_vois(ImageRGBu8::PixelType couleur1, ImageRGBu8::PixelType couleur2){
        couleur1_ = couleur1;
        couleur2_ = couleur2;
        compteur_ = 1;
    }

    color_vois(ImageRGBu8::PixelType couleur1, ImageRGBu8::PixelType couleur2, int i){
        couleur1_ = couleur1;
        couleur2_ = couleur2;
        compteur_ = i;
    }

    color_vois(ImageRGBu8::PixelType couleur1, ImageRGBu8::PixelType couleur2, double i){
        couleur1_ = couleur1;
        couleur2_ = couleur2;
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

    bool operator==(const color_vois& couleur2) const
    {
        return (couleur1_ == couleur2.couleur1_) && (couleur2_ == couleur2.couleur2_);
    }

    bool operator<(const color_vois& couleur2) const
    {
        if(couleur1_ == couleur2.couleur1_){
            return (couleur2_ < couleur2.couleur2_);
        }
        else{
            return (couleur1_ < couleur2.couleur1_);
        }
//        return (couleur1_ < couleur2.couleur1_);// and (couleur2_ < couleur2.couleur2_);
    }
};

bool is_in(std::vector<color_vois>& vecteur, ImageRGBu8::PixelType& couleur1, ImageRGBu8::PixelType& couleur2, int& id)
{
    id = -1;
//    for(ImageRGB8::PixelType col : vecteur)

   const color_vois cv{ couleur1, couleur2 };

    for(int i=0; i<vecteur.size(); i++)
    {
        /*if(vecteur.at(i) == color_vois{couleur1, couleur2})*/
        if (vecteur[i] == cv)
        {
            id = i;
            return true;
        }
    }
    return false;
}


bool is_in(std::vector<color_vois>& vecteur, const color_vois& couleurs, int& id)
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

bool get_doublet(color_vois& doublet, ImageRGBu8::PixelType& couleur1, ImageRGBu8::PixelType& couleur2){
    std::vector couleurs{couleur1, couleur2};
    std::sort(couleurs.begin(), couleurs.end());

    auto last = std::unique(couleurs.begin(), couleurs.end());
    couleurs.erase(last, couleurs.end());

    if(couleurs.size()!=2){return false;}
    else{
        doublet = color_vois{couleurs.at(0), couleurs.at(1)};
        return true;
    }
}


//void update_Tcontent(std::vector<color_vois>& couleurs,
//                     ImageRGBu8::PixelType P1, ImageRGBu8::PixelType P2,
//                     bool i_isEqual_j)
//{
//    int id;
//    bool presence;
//
////    if(P1 != P2){
////        presence = is_in(couleurs, P1, P2, id);
////        if(not presence)
////        {
////            presence = is_in(couleurs, P2, P1, id);
////            if(not presence)
////            {
////                ImageRGBu8::PixelType Pi = std::min(P1,P2);
////                ImageRGBu8::PixelType Pj = std::max(P1,P2);
////                couleurs.push_back(color_vois{Pi, Pj});
////            }
////            else
////            {
////                couleurs.at(id).incr();
////            }
////        }
////        else
////        {
////            couleurs.at(id).incr();
////        }
////    }
//
//
//    if(P1 != P2){
//        presence = is_in(couleurs, P1, P2, id);
//        if(not presence)
//        {
//            couleurs.push_back(color_vois{P1, P2});
//        }
//        else
//        {
//            couleurs.at(id).incr();
//        }
//
//        presence = is_in(couleurs, P2, P1, id);
//        if(not presence)
//        {
//            couleurs.push_back(color_vois{P2, P1});
//        }
//        else
//        {
//            couleurs.at(id).incr();
//        }
//    }
//
//    else if(P1 == P2 and i_isEqual_j){
//        presence = is_in(couleurs, P1, P2, id);
//        if(not presence)
//        {
//            couleurs.push_back(color_vois{P1, P2});
//        }
//        else
//        {
//            couleurs.at(id).incr();
//        }
//    }
//
//
//}


void update_Tcontent(std::vector<color_vois>& couleurs, color_vois& doublet){
    int id = -1;
    bool presence = is_in(couleurs, doublet, id);

    if(! presence){
        couleurs.push_back(doublet);
    }
    else{
        couleurs.at(id).incr();
    }
}

void update_Tcontent(std::vector<color_vois>& couleurs, color_vois& doublet, int qtt){
    int id = -1;
    bool presence = is_in(couleurs, doublet, id);

    if(! presence){
        doublet.count(qtt);
        couleurs.push_back(doublet);
    }
    else{
        couleurs.at(id).incr(qtt);
    }
}

void update_Tcontent(std::vector<color_vois>& couleurs, color_vois& doublet, double qtt){
    int id = -1;
    bool presence = is_in(couleurs, doublet, id);

    if (! presence) {
        doublet.count(qtt);
        couleurs.push_back(doublet);
    }
    else{
        couleurs.at(id).incr(qtt);
    }
}




void update_Hcontent_reel(std::vector<color_vois>& couleurs, ImageRGBu8::PixelType P1, ImageRGBu8::PixelType P2,
                          ImageGrayd histo_, bool i_isEqual_j,
                          int x, int y,  int dx, int dy) {

    int id;
    bool presence;
    double qtt;// = gauss(moy1, moy2, var1, var2, X, Y);
    int cm_size = histo_.height();


    // histo réel
    if(histo_.pixelAbsolute(x, y) < (1./double(cm_size - 1)) || histo_.pixelAbsolute(x + dx, y + dy) < (1./double(cm_size - 1)))
    {
        qtt = 0.;
    }
    else
    {
        qtt = 0.5*(histo_.pixelAbsolute(x, y) + histo_.pixelAbsolute(x + dx, y + dy));
    }


    if(P1 != P2){
        presence = is_in(couleurs, P1, P2, id);
        if(! presence)
        {
            presence = is_in(couleurs, P2, P1, id);
            if (! presence)
            {
                ImageRGBu8::PixelType Pi = std::min(P1,P2);
                ImageRGBu8::PixelType Pj = std::max(P1,P2);
                couleurs.push_back(color_vois{Pi, Pj, qtt});
            }
            else
            {
                couleurs.at(id).incr(qtt);
            }
        }
        else
        {
            couleurs.at(id).incr(qtt);
        }
    }


    else if (P1 == P2 && i_isEqual_j) {

        presence = is_in(couleurs, P1, P2, id);
        if (! presence) {
            couleurs.push_back(color_vois{P1, P2, qtt});
        } else {
            couleurs.at(id).incr(qtt);
        }
    }
}



void update_Hcontent_theo(std::vector<color_vois>& couleurs, ImageRGBu8::PixelType P1, ImageRGBu8::PixelType P2,
                          int cm_size, bool i_isEqual_j,
                          double moy1, double moy2, double var1, double var2,
                          int x, int y,  int dx, int dy) {

    int id;
    bool presence;
    double X;// = (x+0.5)/double(cm_size-1); // évaluation au centre des pixels
    double Y;// = (y+0.5)/double(cm_size-1);
    double qtt;// = gauss(moy1, moy2, var1, var2, X, Y);


    // histo théorique
    X = (x + 0.5 + double(dx) / 2.) / double(cm_size - 1); // évaluation au centre des bord des pixels
    Y = (y + 0.5 + double(dy) / 2.) / double(cm_size - 1);
    qtt = gauss(moy1, moy2, var1, var2, X, Y, double(cm_size - 1));


    if(P1 != P2){
        presence = is_in(couleurs, P1, P2, id);
        if(! presence)
        {
            presence = is_in(couleurs, P2, P1, id);
            if (! presence)
            {
                ImageRGBu8::PixelType P_i = std::min(P1,P2);
                ImageRGBu8::PixelType P_j = std::max(P1,P2);
                couleurs.push_back(color_vois{P_i, P_j, qtt});
            }
            else
            {
                couleurs.at(id).incr(qtt);
            }
        }
        else
        {
            couleurs.at(id).incr(qtt);
        }
    }


    else if (P1 == P2 && i_isEqual_j) {

        presence = is_in(couleurs, P1, P2, id);
        if (! presence) {
            couleurs.push_back(color_vois{P1, P2, qtt});
        } else {
            couleurs.at(id).incr(qtt);
        }
    }
}


#endif //ASTEX_MESURE_VOISINAGE_H
