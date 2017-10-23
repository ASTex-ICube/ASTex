// Module Fragmentation by Kenneth

#ifndef FRAGMENTER_H_
#define FRAGMENTER_H_

#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <limits>
#include <cmath>
#include <assert.h>
#include <string>
#include <sstream>

#include "hvPicture.h"

namespace Fragmentation
{
using namespace hview;

enum TYPE_OF_IMAGE
{
	REGION_INDEX = 1, // bit 1
	BORDERS = 2, // bit 2
	PASTEL_COLORS = 4, // bit 3
	NON_SEGMENTED = 8 // bit 4
} ;

/*
 *	\class Fragmenter managing image fragmentation given three main parameters:
 *	min, max fragment size and threshold value defining when pixels should be separated into two fragments.
 */
class Fragmenter
{
private:
	enum pixelstatus
	{
		NOT_VISITED = 0,
		VISITED  = -1,
		BORDER  = -2,
	} ;

	typedef struct
	{
		unsigned int red, green, blue ;
	} Color ;

	class Matrix
	{
	private:
		const unsigned int m_h, m_w ;
		int *m_data ;

	public:
		Matrix(unsigned int h, unsigned int w):
			m_h(h),
			m_w(w)
		{
			m_data = new int[m_h*m_w] ;
			for (int i=0; i<m_h*m_w; i++) m_data[i]=0;
		}
		~Matrix()
		{
			delete[] m_data ;
		}

		int& operator()(unsigned int h, unsigned int w)
		{
			return m_data[m_w*h + w] ;
		}

		const int& operator()(unsigned int h, unsigned int w) const
		{
			return m_data[m_w*h + w] ;
		}
	};

	bool m_isFragmented ;
	bool m_isTorus ;

	const hvPictRGB<unsigned char> *m_img ;
	Matrix m_status ;
	std::vector<std::vector<std::pair<int,int> > > m_fragments ;

public:
	/**
	 * \brief constructor stores the given image and initializes internal parameters
	 *
	 * \param a valid QImage
	 */
	Fragmenter(const hvPictRGB<unsigned char> *image, bool torus=false) :
			m_isTorus(torus),
			m_isFragmented(false),
			m_img(image),
			m_status(Matrix(m_img->sizeX(),m_img->sizeY()))
	{
		//srand(time(NULL)) ;
	}

	~Fragmenter()
	{}

	int nFragments() const
	{
		if (!m_isFragmented)
		{
			hvFatal("Fragmenter::nothing has been fragmented yet. Did you forget to call Fragmenter::fragment() ?");
			return -1 ;
		}
		return m_fragments.size();
	}

	int idFragment(int i, int j) const
	{
		if (!m_isFragmented)
		{
			hvFatal("Fragmenter::nothing has been fragmented yet. Did you forget to call Fragmenter::fragment() ?");
			return -1 ;
		}
		return m_status(i,j);
	}


	int fragmentSize(int fid)
	{
		if (!m_isFragmented)
		{
			hvFatal("Fragmenter::Cleanup cannot clean anything up since nothing has been fragmented yet. Did you forget to call Fragmenter::fragment() ?");
			return -1 ;
		}
		return m_fragments[fid].size();
	}

	void pixelFragment(int fid, int idpix, int &x, int &y)
	{
		std::pair<int,int> pp = m_fragments[fid][idpix];
		x=pp.first; y=pp.second;
	}
	/**
	 * \brief cuts the given image into fragments by a region-growing algorithm.
	 *
	 * \param max the maximal fragment size : up to 10000
	 * \param threshold the color distance over which neighboring pixels will be separated: between 30 and 100
	 *
	 * \return the amount of fragments
	 */
unsigned int fragment(unsigned int max, unsigned int threshold, hvPict<unsigned char> *pclasses=0)
{
	m_fragments.clear() ;

	std::list<std::pair<int,int> > border ; // collection de pixels de bordure de régions (qui auront été exclus d'une région)
	std::vector<std::pair<int,int> > cur_region ; // collection de pixels de la région en cours de construction
	std::vector<std::pair<int,int> > cur_region_candidates ; // collection de candidats à la région en cours de construction

	// prendre un pixel au hasard
	std::pair<int,int> index(rand() % m_img->sizeX(), rand() % m_img->sizeY()) ;
	// et avancer jusqu'à ce qu'on tombe sur un pixel non encore traité (utile que dans certains cas)
	do
	{
		index.first += 1 ;
		if (index.first >= m_img->sizeX())
		{
			index.first = 0 ;
			++index.second ;
			if (index.second >= m_img->sizeY())
				break ;
		}
	} while(m_status(index.first, index.second) != NOT_VISITED) ;
	if (m_status(index.first, index.second) != NOT_VISITED) hvFatal("could not find pixel!");

	// ajouter ce pixel aux candidats à la région courante
	cur_region_candidates.push_back(index) ;
	int cur_reg_num = 1 ; // init numéro de région

	Color sum = {0, 0, 0} ; // variable qui va contenir la couleur moyenne de la zone courante

//	printf("Fragmenting: starting with %d,%d\n", index.first, index.second);
	// Boucle globale
	while(!m_isFragmented)
	{
		if (cur_region_candidates.empty())
		{
			if (border.empty()) // on s'arrete lorsque tous les pixels auront été traités: plus de candidats à la région courante + plus de pixels qui auront été précédemment exclus d'une région
				m_isFragmented = true ;
			else
			{	// on entame une nouvelle région s'il n'y a plus de candidats à la région courante
				cur_region_candidates.push_back(border.front()) ; // on ajoute le premier parmi les pixels précédemment exclus d'une autre zone
				border.pop_front() ;
				sum.red = 0 ; sum.green = 0 ; sum.blue = 0 ; // re-init de la moyenne
				++cur_reg_num ; // incr. numéro de la zone courante
				//std::cout << "New regions " << cur_reg_num << "of size " << cur_region.size() << std::endl ;
			}
			m_fragments.push_back(cur_region) ; // ajout de la région terminée à la liste globale des régions
			cur_region.clear() ; // init de la région courante
		}
		else
		{	// sinon évaluer les candidats à la région courante
			//std::cout << cur_region_candidates.size() << std::endl ;
			//std::cout << cur_region.size() << std::endl ;
			//std::cout << cur_region_candidates.front().first << ", " << cur_region_candidates.front().second << std::endl ;

			int index = rand()%(cur_region_candidates.size()) ;			// prendre un candidat au hasard
			const std::pair<int,int> current = cur_region_candidates[index] ;

			m_status(current.first, current.second) = cur_reg_num ;	// l'associer à la région courante
			cur_region.push_back(current) ;

			hvColRGB<unsigned char> col = m_img->get(current.first,current.second) ;	// compatbiliser ce pixel dans la moyenne de la région courante
			sum.red += col.RED() ;
			sum.green += col.GREEN() ;
			sum.blue += col.BLUE() ;

			unsigned char label=0; 
			if (pclasses!=0) label = pclasses->get(current.first,current.second) ;

			// évaluer ses quatre voisins et les ajouter soit au bord exclu, soit aux candidats
			int xx,yy;
			//printf("candidate %d(%d)/%d: %d,%d (%d/%d)\n", cur_region.size(), border.size(), max,current.first, current.second,index,cur_region_candidates.size());
			xx=current.first - 1; if (m_isTorus && xx < 0) xx += m_img->sizeX(); 
			std::pair<int,int> next(xx,current.second) ;
			//printf("check ."); fflush(stdout);
			check_neighbor(next, sum, threshold, border, cur_region_candidates, cur_reg_num, cur_region.size(), pclasses, label) ;
			xx=current.first + 1; if (m_isTorus && xx>= m_img->sizeX()) xx -= m_img->sizeX();
			next = std::pair<int,int>(xx,current.second) ;
			//printf("."); fflush(stdout);
			check_neighbor(next, sum, threshold, border, cur_region_candidates, cur_reg_num, cur_region.size(), pclasses, label) ;
			//printf("*(%d,%d)",index,cur_region_candidates.size()); fflush(stdout);
			yy=current.second - 1; if (m_isTorus && yy < 0) yy += m_img->sizeY(); 
			//printf("="); fflush(stdout);
			next = std::pair<int,int>(current.first,yy) ;
			//printf("."); fflush(stdout);
			check_neighbor(next, sum, threshold, border, cur_region_candidates, cur_reg_num, cur_region.size(), pclasses, label) ;
			yy=current.second + 1; if (m_isTorus && yy>= m_img->sizeY()) yy -= m_img->sizeY();
			next = std::pair<int,int>(current.first,yy) ;
			//printf("."); fflush(stdout);
			check_neighbor(next, sum, threshold, border, cur_region_candidates, cur_reg_num, cur_region.size(), pclasses, label) ;
			//printf("\n");
			cur_region_candidates.erase(cur_region_candidates.begin() + index) ; // supprimer le pixel qui vient d'etre ajouté de la liste des candidats

			if (cur_region.size() == max) // si region courante atteint tmax
			{
				while (!cur_region_candidates.empty())
				{
					// migrer les pixels de la liste des candidats à la liste des exclus au bord
					border.push_back(cur_region_candidates.back()) ;
					cur_region_candidates.pop_back() ;
				}
			}
		}
	}

	return m_fragments.size() ;
}
	/**
	 * \brief Removes zones that are to small by dispatching its pixels into another zone.
	 *
	 * Warning : no check is done on previously defined max size used for fragmentation.
	 *
	 * \param min the minimal value
	 *
	 * \return the amount of resulting fragments
	 */
unsigned int cleanup(unsigned int min)
{
	if (!m_isFragmented)
	{
		hvFatal("Fragmenter::Cleanup cannot clean anything up since nothing has been fragmented yet. Did you forget to call Fragmenter::fragment() ?");
		return -1 ;
	}

	// crée une map ordonnée de regions avec leur taille
	std::map<std::pair<unsigned int,int>, std::vector<std::pair<int,int> >* > region_order ; // triplet (taille, num frag, frag)
	for (unsigned int i = 0 ; i < m_fragments.size() ; ++i )
	{
		// insert each region with associated size and region number
		std::pair<std::pair<unsigned int,int>, std::vector<std::pair<int,int> > *> reg_i ;
		reg_i.first = std::pair<unsigned int,int>((unsigned int)(m_fragments[i].size()), m_status(m_fragments[i][0].first, m_fragments[i][0].second)) ; // couple (taille, num frag)
		reg_i.second = &m_fragments[i] ; // (pointeur sur fragment)
		region_order.insert(reg_i) ;
	}

	for (std::map<std::pair<unsigned int,int>, std::vector<std::pair<int,int> >* >::iterator it = region_order.begin() ;	// iterateur sur la map
			it->first.first < min && it != region_order.end() ; // condition d'arrêt : aucune région de taille trop petite
			it = region_order.begin()) // remise sur le premier élément de la liste à chaque itération
	{
		std::vector<std::pair<int,int> > *cur_reg = it->second ; // get la région trop petite à vider

		// ajouter tout ses pixels dans une liste
		std::list<std::pair<int,int> > pixels ;
		for (std::vector<std::pair<int,int> >::const_iterator pix = cur_reg->begin() ; pix != cur_reg->end() ; ++pix)
		{
			pixels.push_back(*pix) ;
		}

		// vide la liste de pixels en les affectant à d'autres zones
		while (!pixels.empty())
		{
			// créer la liste des régions voisines candidates à la reprise
			std::vector<std::pair<int,int> > candidates ;
			int xx,yy;
			xx=pixels.front().first - 1; if (m_isTorus && xx < 0) xx += m_img->sizeX(); 
			std::pair<int,int> neighbor ;
			neighbor = std::pair<int,int>(xx, pixels.front().second) ;
			if (xx >= 0 && m_status(pixels.front().first,pixels.front().second) != m_status(neighbor.first,neighbor.second))
				candidates.push_back(neighbor) ;
			xx=pixels.front().first + 1; if (m_isTorus && xx >= m_img->sizeX()) xx -= m_img->sizeX(); 
			neighbor = std::pair<int,int>(xx, pixels.front().second) ;
			if (xx < m_img->sizeX()  && m_status(pixels.front().first,pixels.front().second) != m_status(neighbor.first,neighbor.second))
				candidates.push_back(neighbor) ;
			yy=pixels.front().second - 1; if (m_isTorus && yy < 0) yy += m_img->sizeY(); 
			neighbor = std::pair<int,int>(pixels.front().first, yy) ;
			if (yy >= 0 && m_status(pixels.front().first,pixels.front().second) != m_status(neighbor.first,neighbor.second))
				candidates.push_back(neighbor) ;
			yy=pixels.front().second + 1; if (m_isTorus && yy >= m_img->sizeY()) yy -= m_img->sizeY(); 
			neighbor = std::pair<int,int>(pixels.front().first, yy) ;
			if (yy < m_img->sizeY()  && m_status(pixels.front().first,pixels.front().second) != m_status(neighbor.first,neighbor.second))
				candidates.push_back(neighbor) ;

			// si aucun pixel voisin n'est d'une autre région
			if (candidates.empty())
				pixels.push_back(pixels.front()) ; // remettre le pixel en fond de list
			else
			{
				 //std::cout << pixels.size() << std::endl ;

				// Trouver le fragment le moins éloigné parmi les candidats
				Color px ;
				hvColRGB<unsigned char> cc=m_img->get(pixels.front().first, pixels.front().second);
				px.red = cc.RED() ;
				px.green = cc.GREEN() ;
				px.blue = cc.BLUE() ;

				int minid = 0 ;
				float min = 999 ;
				for (unsigned int i = 0 ; i < candidates.size() ; ++i)
				{
					Color candidate ;
					hvColRGB<unsigned char> pcc=m_img->get(candidates[i].first, candidates[i].second);
					candidate.red = pcc.RED() ;
					candidate.green = pcc.GREEN() ;
					candidate.blue = pcc.BLUE() ;
					candidate.red -= px.red ;	candidate.green -= px.green ;	candidate.blue -= px.blue ;
					float diff = colorNorm(candidate) ;
					if (diff < min)
					{
						min = diff ; // distance min entre pixel voisin et pixel courant
						minid = i ; // id du fragment le moins éloigné
					}
				}
				int no_neigh = m_status(candidates[minid].first, candidates[minid].second) ; // numéro de zone du meilleur voisin
				unsigned int old_size = m_fragments[no_neigh - 1].size() ; // taille fragment du meilleur voisin

				// get region voisine
				std::map<std::pair<unsigned int,int>, std::vector<std::pair<int,int> >* >::iterator it_neigh ;
				it_neigh = region_order.find(std::pair<unsigned int,int>(old_size, no_neigh)) ;
				if (it_neigh == region_order.end())
					hvFatal("something went wrong!");

				// ajout du pixel courant à la zone voisinne
				it_neigh->second->push_back(pixels.front()) ;
				m_status(pixels.front().first,pixels.front().second) = no_neigh ;

				// ajout voisin mis à jour ds map des regions
				std::pair<std::pair<unsigned int,int>, std::vector<std::pair<int,int> > *> newreg_neigh ;
				newreg_neigh.second = it_neigh->second ;
				newreg_neigh.first.first = it_neigh->second->size() ;
				newreg_neigh.first.second = no_neigh ;
				region_order.insert(newreg_neigh) ;

				// supprime ancien voisin de la map des regions
				region_order.erase(std::pair<int,int>(old_size,no_neigh)) ;
			}
			pixels.pop_front() ; // supprimer le pixel, qu'il ait été traité ou remis en fond de liste
		}

		 // supprime la region vidée de la map
		region_order.erase(it) ;
	}

	// reset m_fragments (non maintenu lors du nettoyage)
	std::vector<std::vector<std::pair<int,int> > > fragments ;
	fragments.reserve(region_order.size()) ;
	for (std::map<std::pair<unsigned int,int>, std::vector<std::pair<int,int> >* >::const_iterator it = region_order.begin() ;	// iterateur sur la map
				it != region_order.end() ;
				it++)
	{
		fragments.push_back(*it->second) ;
	}
	m_fragments = fragments ;

	return m_fragments.size() ;
}

	/**
	 * \brief saves the image to a file.
	 *
	 * \param name the path and name of the file (will be appended with type and extension)
	 * \param type one or more of the TYPE_OF_IMAGE values (can be summed with a logical OR).
	 * Usage : saveImageToPNG("/home/myfile", REGION_INDEX | BORDER) const ;
	 */

void toPictRGB(const hvPict<unsigned char> *featid, hvPictRGB<unsigned char> &res, unsigned int type) const
{
	// compute highest fragment number
	int no_max = 0 ;
	for( int i=0; i<m_img->sizeX(); ++i )
		for( int j=0; j<m_img->sizeY(); ++j )
			if( m_status(i,j) > no_max )
				no_max = m_status(i,j);

	if (no_max > 256*256)
		hvFatal("Warning: there are more fragments than possible colors");

	if ((type & REGION_INDEX) && m_isFragmented)
	{
		// create array of unique colors
		//hvColRGB<unsigned char> *colors = new hvColRGB<unsigned char>[no_max+1] ;
		//for ( int i = 0 ; i <= no_max ; ++i)
		//{
			//colors[i]=hvColRGB<unsigned char>((unsigned char)(i%256), (unsigned char)((i / 256) % 256), (unsigned char)((i / (256*256)) % 256));
			//colors[i].setRed  ( i % 256 ) ;
			//colors[i].setGreen( (i / 256) % 256 ) ;
			//colors[i].setBlue ( (i / (256*256)) % 256 ) ;
		//}

		// export an image where each fragment has its own unique colors.
		res.reset(m_img->sizeX(), m_img->sizeY(),  hvColRGB<unsigned char>(0)) ;
		for (int i = 0 ; i < m_img->sizeX() ; ++i)
		{
			for (int j = 0 ; j < m_img->sizeY() ; ++j)
			{
				int id=m_status(i,j);
				int ii=(featid==0?0:(int)featid->get(i,j)+1)*53;
				hvColRGB<unsigned char> col((unsigned char)((id*7+ii)%128+128),(unsigned char)(((id*11+ii) / 256) % 128+128),(unsigned char)(ii%128+128));
				res.update(i,j,col) ;
			}
		}

		//delete[] colors ;
	}
}

void toPictRGB(hvPictRGB<unsigned char> &res) const
{
	// compute highest fragment number
	int no_max = 0 ;
	for( int i=0; i<m_img->sizeX(); ++i )
		for( int j=0; j<m_img->sizeY(); ++j )
			if( m_status(i,j) > no_max )
				no_max = m_status(i,j);

	if (no_max > 256*256*256)
		hvFatal("Warning: there are more fragments than possible colors");

	//printf("highest color: %d\n", no_max);

	if (m_isFragmented)
	{
		// create array of unique colors
		hvColRGB<unsigned char> *colors = new hvColRGB<unsigned char>[no_max+1] ;
		for ( int i = 0 ; i <= no_max ; ++i)
		{
			colors[i]=hvColRGB<unsigned char>((unsigned char)(i%256), (unsigned char)((i / 256) % 256), (unsigned char)((i / (256*256)) % 256));
			//colors[i].setRed  ( i % 256 ) ;
			//colors[i].setGreen( (i / 256) % 256 ) ;
			//colors[i].setBlue ( (i / (256*256)) % 256 ) ;
		}

		// export an image where each fragment has its own unique colors.
		res.reset(m_img->sizeX(), m_img->sizeY(),  hvColRGB<unsigned char>(0)) ;
		for (int i = 0 ; i < m_img->sizeX() ; ++i)
		{
			for (int j = 0 ; j < m_img->sizeY() ; ++j)
			{
				int id=m_status(i,j);
				//int ii=(featid==0?0:(int)featid->get(i,j)+1)*53;
				//hvColRGB<unsigned char> col((unsigned char)((id*7+ii)%128+128),(unsigned char)(((id*11+ii) / 256) % 128+128),(unsigned char)(ii%128+128));
				res.update(i,j,colors[id]) ;
			}
		}

		delete[] colors ;
	}
}
	/**
	 * \brief returns the amount of fragments currently stored
	 *
	 * \return the amount of fragments
	 */
	unsigned int getNbFragments() const
	{
		return m_fragments.size() ;
	}

void check_neighbor(const std::pair<int,int>& next, const Color& sum, const float& tol, std::list<std::pair<int,int> >& border, 
std::vector<std::pair<int,int> >& cur_region_candidates, int cur_reg_num, int cur_reg_elems, hvPict<unsigned char> *pclasses=0, unsigned char label=0)
{
	// vérif si le pixel est dans l'image
	if ((next.first < 0 || next.second < 0 || next.first >= m_img->sizeX() || next.second >= m_img->sizeY()))
		return ;

	// vérif si le pixel a bien jamais été visité
	if (m_status(next.first,next.second) != NOT_VISITED)
		return ;

	// get couleur du pixel à évaluer
	hvColRGB<unsigned char> col = m_img->get(next.first,next.second) ;
	Color myColor = { (unsigned int)col.RED() - sum.red / cur_reg_elems, (unsigned int)col.GREEN() - sum.green / cur_reg_elems, (unsigned int)col.BLUE() - sum.blue / cur_reg_elems }  ;

	unsigned char nlabel=0;
	if (pclasses!=0) nlabel = pclasses->get(next.first,next.second) ;

	// vérifie si ce pixel excède le seuil de tolérance p/r à la moyenne
	if (nlabel==label && colorNorm(myColor) < tol)
	{
		// ajout
		m_status(next.first, next.second) = cur_reg_num ;
		cur_region_candidates.push_back(next) ;
		//std::cout << "region add" << std::endl ;
	}
	else
	{
		// exclusion (= ajout au bord)
		m_status(next.first, next.second) = BORDER ;
		border.push_back(next) ;
		//std::cout << "border add" << std::endl ;
	}
}

	static float colorNorm(const Color& c)
	{
		return sqrt((float)(c.red * c.red + c.green * c.green + c.blue * c.blue)) ;
	}


} ;

} // namespace Fragmentation

#endif /* FRAGMENTER_H_ */
