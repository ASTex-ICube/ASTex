
#include "Fragmenter.h"

namespace Fragmentation
{

Fragmenter::Fragmenter(const QImage& image) :
		m_isFragmented(false),
		m_img(image),
		m_status(Matrix(m_img.height(),m_img.width()))
{
	srand(time(NULL)) ;
}

Fragmenter::~Fragmenter()
{}

unsigned int
Fragmenter::fragment(unsigned int max, unsigned int threshold)
{
	m_fragments.clear() ;

	std::list<std::pair<int,int> > border ; // collection de pixels de bordure de régions (qui auront été exclus d'une région)
	std::vector<std::pair<int,int> > cur_region ; // collection de pixels de la région en cours de construction
	std::vector<std::pair<int,int> > cur_region_candidates ; // collection de candidats à la région en cours de construction

	// prendre un pixel au hasard
	std::pair<int,int> index(rand() % m_img.height(), rand() % m_img.width()) ;
	// et avancer jusqu'à ce qu'on tombe sur un pixel non encore traité (utile que dans certains cas)
	do
	{
		index.first += 1 ;
		if (index.first == m_img.height())
		{
			index.first = 0 ;
			++index.second ;
			if (index.second == m_img.width())
				break ;
		}
	} while(m_status(index.first, index.second) != NOT_VISITED) ;

	// ajouter ce pixel aux candidats à la région courante
	cur_region_candidates.push_back(index) ;
	int cur_reg_num = 1 ; // init numéro de région

	Color sum = {0, 0, 0} ; // variable qui va contenir la couleur moyenne de la zone courante

	// Boucle globale
	while(!m_isFragmented)
	{
		if (cur_region_candidates.empty())
		{
			if (border.empty()) // on s'arrete lorsque tous les pixels auront été traités: plus de candidats à la région courante + plus de pixels qui auront été précédemment exclus d'une région
				m_isFragmented = true ;
			else
			{	// on entame une nouvelle région s'il n'y a plus de candidats à la région courante
				//std::cout << "New regions of size " << cur_region.size() << std::endl ;
				cur_region_candidates.push_back(border.front()) ; // on ajoute le premier parmi les pixels précédemment exclus d'une autre zone
				border.pop_front() ;
				sum.red = 0 ; sum.green = 0 ; sum.blue = 0 ; // re-init de la moyenne
				++cur_reg_num ; // incr. numéro de la zone courante
			}
			m_fragments.push_back(cur_region) ; // ajout de la région terminée à la liste globale des régions
			cur_region.clear() ; // init de la région courante
		}
		else
		{	// sinon évaluer les candidats à la région courante
//			std::cout << cur_region_candidates.size() << std::endl ;
//			std::cout << cur_region.size() << std::endl ;
//			std::cout << cur_region_candidates.front().first << ", " << cur_region_candidates.front().second << std::endl ;

			int index = rand()%(cur_region_candidates.size()) ;			// prendre un candidat au hasard
			const std::pair<int,int>& current = cur_region_candidates[index] ;

			m_status(current.first, current.second) = cur_reg_num ;	// l'associer à la région courante
			cur_region.push_back(current) ;

			QRgb col = m_img.pixel(current.first,current.second) ;	// compatbiliser ce pixel dans la moyenne de la région courante
			sum.red += qRed(col) ;
			sum.green += qGreen(col) ;
			sum.blue += qBlue(col) ;

			// évaluer ses quatre voisins et les ajouter soit au bord exclu, soit aux candidats
			std::pair<int,int> next(current.first - 1,current.second) ;
			check_neighbor(next, sum, threshold, border, cur_region_candidates, cur_reg_num, cur_region.size()) ;
			next = std::pair<int,int>(current.first + 1,current.second) ;
			check_neighbor(next, sum, threshold, border, cur_region_candidates, cur_reg_num, cur_region.size()) ;
			next = std::pair<int,int>(current.first,current.second - 1) ;
			check_neighbor(next, sum, threshold, border, cur_region_candidates, cur_reg_num, cur_region.size()) ;
			next = std::pair<int,int>(current.first,current.second + 1) ;
			check_neighbor(next, sum, threshold, border, cur_region_candidates, cur_reg_num, cur_region.size()) ;

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

unsigned int
Fragmenter::cleanup(unsigned int min)
{
	if (!m_isFragmented)
	{
		std::cerr << "Fragmenter::Cleanup cannot clean anything up since nothing has been fragmented yet. Did you forget to call Fragmenter::fragment() ?" << std::endl ;
		return -1 ;
	}

	// crée une map ordonnée de regions avec leur taille
	std::map<std::pair<unsigned int,int>, std::vector<std::pair<int,int> >* > region_order ; // triplet (taille, num frag, frag)
	for (unsigned int i = 0 ; i < m_fragments.size() ; ++i )
	{
		// insert each region with associated size and region number
		std::pair<std::pair<unsigned int,int>, std::vector<std::pair<int,int> > *> reg_i ;
		reg_i.first = std::pair<unsigned int,int>(m_fragments[i].size(), m_status(m_fragments[i][0].first, m_fragments[i][0].second)) ; // couple (taille, num frag)
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
			std::pair<int,int> neighbor ;
			neighbor = std::pair<int,int>(pixels.front().first - 1, pixels.front().second) ;
			if (pixels.front().first > 0 && m_status(pixels.front().first,pixels.front().second) != m_status(neighbor.first,neighbor.second))
				candidates.push_back(neighbor) ;
			neighbor = std::pair<int,int>(pixels.front().first + 1, pixels.front().second) ;
			if (pixels.front().first < m_img.height() - 1 && m_status(pixels.front().first,pixels.front().second) != m_status(neighbor.first,neighbor.second))
				candidates.push_back(neighbor) ;
			neighbor = std::pair<int,int>(pixels.front().first, pixels.front().second - 1) ;
			if (pixels.front().second > 0 && m_status(pixels.front().first,pixels.front().second) != m_status(neighbor.first,neighbor.second))
				candidates.push_back(neighbor) ;
			neighbor = std::pair<int,int>(pixels.front().first, pixels.front().second + 1) ;
			if (pixels.front().second < m_img.width() - 1 && m_status(pixels.front().first,pixels.front().second) != m_status(neighbor.first,neighbor.second))
				candidates.push_back(neighbor) ;

			// si aucun pixel voisin n'est d'une autre région
			if (candidates.empty())
				pixels.push_back(pixels.front()) ; // remettre le pixel en fond de list
			else
			{
				 //std::cout << pixels.size() << std::endl ;

				// Trouver le fragment le moins éloigné parmi les candidats
				Color px ;
				px.red = qRed(m_img.pixel(pixels.front().first, pixels.front().second)) ;
				px.green = qGreen(m_img.pixel(pixels.front().first, pixels.front().second)) ;
				px.blue = qBlue(m_img.pixel(pixels.front().first, pixels.front().second)) ;

				int minid = 0 ;
				float min = 999 ;
				for (unsigned int i = 0 ; i < candidates.size() ; ++i)
				{
					Color candidate ;
					candidate.red = qRed(m_img.pixel(candidates[i].first, candidates[i].second)) ;
					candidate.green = qGreen(m_img.pixel(candidates[i].first, candidates[i].second)) ;
					candidate.blue = qBlue(m_img.pixel(candidates[i].first, candidates[i].second)) ;
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
					std::cout << "something went wrong!" << std::endl ;

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

void
Fragmenter::check_neighbor(const std::pair<int,int>& next, const Color& sum, const float& tol, std::list<std::pair<int,int> >& border, std::vector<std::pair<int,int> >& cur_region_candidates, int cur_reg_num, int cur_reg_elems)
{
	// vérif si le pixel est dans l'image
	if (next.first < 0 || next.second < 0 || next.first >= m_img.height() || next.second >= m_img.width())
		return ;

	// vérif si le pixel a bien jamais été visité
	if (m_status(next.first,next.second) != NOT_VISITED)
		return ;

	// get couleur du pixel à évaluer
	QRgb col = m_img.pixel(next.first, next.second) ;
	Color myColor = { qRed(col) - sum.red / cur_reg_elems, qGreen(col) - sum.green / cur_reg_elems, qBlue(col) - sum.blue / cur_reg_elems }  ;

	// vérifie si ce pixel excède le seuil de tolérance p/r à la moyenne
	if (colorNorm(myColor) < tol)
	{
		// ajout
		m_status(next.first, next.second) = cur_reg_num ;
		cur_region_candidates.push_back(next) ;
//		std::cout << "region add" << std::endl ;
	}
	else
	{
		// exclusion (= ajout au bord)
		m_status(next.first, next.second) = BORDER ;
		border.push_back(next) ;
//		std::cout << "border add" << std::endl ;
	}
}

void
Fragmenter::saveImageToPNG(const char* string, unsigned int type) const
{
	// compute highest fragment number
	int no_max = 0 ;
	for( int i=0; i<m_img.height(); ++i )
		for( int j=0; j<m_img.width(); ++j )
			if( m_status(i,j) > no_max )
				no_max = m_status(i,j);

	if (no_max > 256*256*256)
		std::cerr << "Warning: there are more fragments than colors" << std::endl ;

	if ((type & REGION_INDEX) && m_isFragmented)
	{
#ifdef DISPLAY
		std::cout << "Exporting region indices" << std::endl ;
#endif
		// create array of unique colors
		QColor *colors = new QColor[no_max+1] ;
		for ( int i = 0 ; i <= no_max ; ++i)
		{
			colors[i].setRed  ( i % 256 ) ;
			colors[i].setGreen( (i / 256) % 256 ) ;
			colors[i].setBlue ( (i / (256*256)) % 256 ) ;
		}

		// export an image where each fragment has its own unique colors.
		QImage indexImage(m_img.height(), m_img.width(),  m_img.format()) ;
		for (int i = 0 ; i < m_img.height() ; ++i)
		{
			for (int j = 0 ; j < m_img.width() ; ++j)
			{
				indexImage.setPixel(i,j,colors[m_status(i,j)].rgb()) ;
			}
		}

		delete[] colors ;

		std::stringstream s ;
		s << string << "_index.png" ;
		indexImage.save(s.str().c_str(), NULL, 100) ;
	}
	if ((type & BORDERS) && m_isFragmented)
	{
#ifdef DISPLAY
		std::cout << "Exporting fragment borders" << std::endl ;
#endif
		// Create a binary image where the fragment borders are black.
		QImage regionsImage(m_img.height(), m_img.width(),  m_img.format()) ;
		for (int i = 0 ; i < m_img.height() ; ++i)
		{
			for (int j = 0 ; j < m_img.width() ; ++j)
			{
				regionsImage.setPixel(i,j,qRgb(255,255,255)) ;
				if (m_status(i,j) != m_status(std::max(i-1,0),j))
					regionsImage.setPixel(i,j,qRgb(0,0,0)) ;
				if (m_status(i,j) != m_status(i,std::max(j-1,0)))
					regionsImage.setPixel(i,j,qRgb(0,0,0)) ;
				if (m_status(i,j) != m_status(std::min(i+1,m_img.height()-1),j))
					regionsImage.setPixel(i,j,qRgb(0,0,0)) ;
				if (m_status(i,j) != m_status(i,std::min(j+1,m_img.width()-1)))
					regionsImage.setPixel(i,j,qRgb(0,0,0)) ;
			}
		}

		std::stringstream s ;
		s << string << "_region.png" ;
		regionsImage.save(s.str().c_str(), NULL, 100) ;
	}
	if ((type & PASTEL_COLORS) && m_isFragmented)
	{
#ifdef DISPLAY
		std::cout << "Exporting pastel colors" << std::endl ;
#endif

		// color each zone with random pastel colors.
		QColor *colors = new QColor[no_max+1] ;
		for ( int i = 0 ; i <= no_max ; ++i)
		{
			colors[i].setRed  ( 128 + (rand()&127) );
			colors[i].setGreen( 128 + (rand()&127) );
			colors[i].setBlue ( 128 + (rand()&127) );
		}

		QImage pastelImage(m_img.height(), m_img.width(),  m_img.format()) ;
		for (int i = 0 ; i < m_img.height() ; ++i)
		{
			for (int j = 0 ; j < m_img.width() ; ++j)
			{
				pastelImage.setPixel(i,j,colors[m_status(i,j)].rgb()) ;
			}
		}

		delete[] colors ;

		std::stringstream s ;
		s << string << "_pastel.png" ;
		pastelImage.save(s.str().c_str(), NULL, 100) ;
	}
	if (type & NON_SEGMENTED)
	{
#ifdef DISPLAY
		std::cout << "Exporting non-segmented colors" << std::endl ;
#endif
		std::stringstream s ;
		s << string << ".png" ;
		m_img.save(s.str().c_str(), NULL, 100) ;
	}
}

} // namespace Fragmenter

