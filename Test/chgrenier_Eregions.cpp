//
// Created by grenier on 10/08/23.
//

#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>

using namespace ASTex;

struct Region_info{
public:
    int id_;
    double couleur_;
    std::vector<std::array<int, 2>> pixel_contenus_;
    std::vector<std::array<int, 2>> pixel_voisins_;

    Region_info(int id, double couleur, std::array<int, 2>& pixel_contenus, std::vector<std::array<int, 2>>& pixel_voisins)
    : id_(id), couleur_(couleur), pixel_voisins_(pixel_voisins)
    {
        pixel_contenus_.push_back(pixel_contenus);
    }

    void add_contenu(std::array<int, 2> pixel){
        pixel_contenus_.push_back(pixel);
    }

    void add_voisin(std::array<int, 2> pixel){
        pixel_voisins_.push_back(pixel);
    }
};

template <class Region_info>
class Graph_regions{
    std::map<int, Region_info> sommets_;
    ImageGrayu8 input_;

    std::vector<std::array<int, 2>> listing_voisins(std::array<int, 2>& pixel)
    {
        std::vector<std::array<int, 2>> voisin;

        if(pixel.at(0) -1 >= 0){
            voisin.push_back(std::array<int, 2>{pixel.at(0) -1, pixel.at(1)});
        }
        if(pixel.at(1) -1 >= 0){
            voisin.push_back(std::array<int, 2>{pixel.at(0), pixel.at(1) -1});
        }

        if(pixel.at(0) +1 < input_.width()){
            voisin.push_back(std::array<int, 2>{pixel.at(0) +1, pixel.at(1)});
        }
        if(pixel.at(1) +1 < input_.height()){
            voisin.push_back(std::array<int, 2>{pixel.at(0), pixel.at(1) +1});
        }

        return voisin;
    }

    bool is_in(std::array<int, 2>& pixel, std::vector<std::array<int, 2>>& pixel_list){
        return std::find(pixel_list.begin(), pixel_list.end(), pixel) != pixel_list.end();
    }

    std::array<int, 2> new_ancrage(){
        // TODO
    }


public:
    Graph_regions(ImageGrayu8& image) : input_(image)
    {

    }

    void create_graph_regions(){
        bool fini = false;
        int id_sommet = 0;
        int x_coord = 0;
        int y_coord = 0;
        std::array<int, 2> suivant{x_coord, y_coord};

        while (not fini){
            bool tous_voisins_teste = false;
            int id_voisin = 0;

            sommets_[id_sommet] = Region_info{id_sommet, input_.pixelAbsolute(x_coord, y_coord), suivant, listing_voisins(suivant)};

            while (not tous_voisins_teste){
                x_coord = sommets_[id_sommet].pixel_voisins_.at(id_voisin).at(0);
                y_coord = sommets_[id_sommet].pixel_voisins_.at(id_voisin).at(1);

                if (sommets_[id_sommet].couleur_ == input_.pixelAbsolute(x_coord, y_coord)){ // voisin testé de la même couleur
                    sommets_[id_sommet].add_contenu(std::array<int, 2>{x_coord, y_coord}); // on ajoute le pixel dans la région
                    sommets_[id_sommet].pixel_voisins_.erase(id_voisin); // on retire le pixel de la liste des voisins

                    std::vector<std::array<int, 2>> nouveaux_voisins = listing_voisins(std::array<int, 2>{x_coord, y_coord});
                    for(auto voisin : nouveaux_voisins){
                        if (!is_in(voisin, sommets_[id_sommet].pixel_voisins_) and !is_in(voisin, sommets_[id_sommet].pixel_contenus_)){
                            sommets_[id_sommet].add_voisin(voisin);
                        }
                    }

                    if(id_voisin >= sommets_[id_sommet].pixel_voisins_.size()){
                        tous_voisins_teste = true;
                    }
                }
                else {
                    id_voisin += 1;
                    if (id_voisin >= sommets_[id_sommet].pixel_voisins_.size()) {
                        tous_voisins_teste = true;
                    }
                }
            }
            id_sommet += 1;
            suivant = new_ancrage();
            if(suivant.at(0) >= input_.width() or suivant.at(1) >= input_.height()){
                fini = true;
            }
        }

    }


};




int main(){
    // ---------------------------------------------------------------------------
    // récupération de l'exemple
    int expl_size = 24;

    ImageGrayu8 exemple_(expl_size, expl_size);
    exemple_.load("/home/grenier/Documents/ASTex_fork/results/T_analysis/data_t.png");


    // graphe des régions
    Graph_regions<Region_info> Regions_exemple{exemple_};
    Regions_exemple.create_graph_regions();


    return 0;
}