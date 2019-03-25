#ifndef __HISTOGRAM__H__
#define __HISTOGRAM__H__

#include <cmath>
#include <ASTex/special_io.h>
#include <ASTex/fourier.h>
#include <ASTex/local_spectrum.h>
#include <ASTex/utils.h>
#include <ASTex/colorspace_filters.h>
#include <ASTex/mask.h>
#include <type_traits>

using namespace ASTex;

/**
 * \brief Histogram base class for building any histogram.
 */
template <class I, class Compare = std::less<typename I::PixelType> >
class Histogram
{
public:

    /**
     * \brief constructs an empty histogram
     */
    Histogram();

    /**
     * \brief constructs a histogram based on the first parameter
     * \param [in] image any pixels container which has members width(), height() a PixelType type and be iterable with for_all_pixels.
     */
    Histogram(const I &image);

    //static

	static bool saveImageToCsv(const I& image, const std::string &out);
	static bool loadImageFromCsv(I& image, const std::string& in);

    //types

    using real = double;
    using PixelType = typename I::PixelType;
	using MapStruct = std::map<PixelType, int, Compare>;

    //iterators

	using iterator = typename MapStruct::iterator;
	using const_iterator = typename MapStruct::const_iterator;

    iterator begin() {return m_histogram.begin();}
    const_iterator begin() const {return m_histogram.begin();}

    iterator end() {return m_histogram.end();}
    const_iterator end() const {return m_histogram.end();}

    iterator find(const PixelType& key) {return m_histogram.find(key);}
    const_iterator find(const PixelType& key) const {return m_histogram.find(key);}

    //others

    /**
     * \brief triggers a routine which updates any statistic stored by the class, such as mean, variance, covariance...
     */
	virtual void updateStatistics();

    /**
     * \brief get the number of pixels in the histogram.
     * \return number of pixels
     */
    int size() const {return m_size;}

    /**
     * \brief get the number of individual bins in the histogram. It is always <= to size().
     * \return number of bins
     */
    int binsNumber() const {return m_histogram.size();}

    /**
     * \brief reset the histogram and builds it with a new image
     * \param [in] image which image the histogram should be built from
     */
    virtual void compute(const I& image);

    /**
     * \brief adds "occ" occurences for the bin "bin" in the histogram
     * \param [in] bin the bin to find or add in the histogram
     * \param [in] occ the number of pixels equal to the bin to be added
     */
    void addPixel(const PixelType &bin, int occ=1);

    /**
     * \brief removes a maximum of "occ" occurences for the bin "bin" in the histogram
     * \param [in] bin the bin to find in the histogram
     * \param [in] occ the number of pixels equal to the bin to be removed
     */
    void removePixel(const PixelType &bin, int occ=1);

    /**
     * \brief sets "occ" occurences for the bin "bin" in the histogram
     * \param [in] bin the bin to find or add in the histogram
     * \param [in] occ the number of pixels equal to the bin to be set
     */
    void setBin(const PixelType &bin, int occ);

    /**
     * \brief completely removes a bin from the histogram
     * \param [in] bin the bin to remove
     */
    void removeBin(const PixelType &bin);

    /**
     * \brief clears the histogram
     */
    virtual void clear();

    /**
     * \brief returns a mean of the pixel type instead of always being a scalar double or array double. Convenient for polymorphism.
     */
    virtual PixelType meanPixelType() const {return PixelType();}

    //out

    /**
     * \brief saves a compact representation of the histogram into a file
     * \param [in] out name of the file to create and write (does not append, be careful)
     */
    virtual void saveHistogram(const std::string& out, int nb_classes_per_dimension=0) const;

    /**
     * \brief saves a sorted list of all pixels in the histogram
     * \param [in] out name of the file to create and write (does not append, be careful)
     */
    virtual void saveFullHistogram(const std::string& out) const;

    //bin to bin comparasions

    /**
     * \brief use histogram intersection to compute a distance between histograms
     * \param [in] other histogram this histogram is compared with
     * \return intersection distance between histograms
     */
    real compareIntersectionWith(const Histogram<I, Compare>& other) const;

    /**
     * \brief use the L-norm Minkowski distance to compute a distance between histograms
     * \param [in] other histogram this histogram is compared with
     * \return L-norm Minkowski distance between histograms
     */
    template <int L_norm=1>
    real compareMinkowskiWith(const Histogram<I, Compare>& other) const;

    /**
     * \brief use the Chi-squared distance to compute a distance between histograms
     * \param [in] other histogram this histogram is compared with
     * \return Chi-squared distance between histograms
     * WARNING : NOT CORRECT
     */
    real compareChi2With(const Histogram<I, Compare>& other) const;

    //tests of fit with distributions

    /**
     * \brief use the Kolmogorov-Smirnov distance to test if this histogram's distribution is uniform
     * \pre PixelType is an arithmetic type (int, float, char...)
     * \param [in] inf expected inferior bound of the histogram and uniform distribution
     * \param [in] sup expected superior bound of the histogram and uniform distribution
     * \return Kolmogorov-Smirnov distance with standard uniform distribution
     */
    template <typename NUMBER=typename I::PixelType, typename std::enable_if<std::is_arithmetic<NUMBER>::value>::type* = nullptr>
    real fitsUniformKS(NUMBER inf=NUMBER(0), NUMBER sup=NUMBER(1)) const;

protected:

    /**
     * \brief Convenient bin to bin comparasion function that uses lambda functions to know what bin to bin test to do
     * \param [in] other histogram this histogram is compared with
     * \param [in] lambdaBothExist (real, real) -> real lambda function to be used when a bin (in frequency) exists in both histograms
     * \param [in] lambdaFirstExists real -> real lambda function to be used when a bin (in frequency) exists only in this histogram
     * \param [in] lambdaSecondExists real -> real lambda function to be used when a bin (in frequency) exists only in the other histogram
     * \return computed bin to bin distance
     */
    template<typename F_bothExist, typename F_firstExists, typename F_secondExists>
    real compareB2BWith(const Histogram<I, Compare>& other,
                        F_bothExist &lambdaBothExist,
                        F_firstExists &lambdaFirstExists,
                        F_secondExists &lambdaSecondExists) const;

	I m_input;
	MapStruct m_histogram;
    int             m_size;
};



//Histogram


template <class I, class Compare>
Histogram<I, Compare>::Histogram() : m_histogram(), m_size(0)
{

}


template <class I, class Compare>
Histogram<I, Compare>::Histogram(const I &image) : m_histogram(), m_size(0), m_input(image)
{
    compute(image);
}

template <class I, class Compare>
void Histogram<I, Compare>::updateStatistics()
{

}

template <class I, class Compare>
bool Histogram<I, Compare>::saveImageToCsv(const I& image, const std::string& out)
{
    std::ofstream ofs_out(out);
    ofs_out << image.width() << std::endl;
    ofs_out << image.height() << std::endl;
    image.for_all_pixels([&] (const typename I::PixelType &pix)
    {
        ofs_out << pix << std::endl;
    });
    ofs_out.close();
	if(!ofs_out)
		return false;
	return true;
}

template <class I, class Compare>
bool Histogram<I, Compare>::loadImageFromCsv(I& image, const std::string &in)
{
	const unsigned pixelSize = sizeof(typename I::PixelType)/sizeof(typename I::DataType);
    std::ifstream ifs_in(in);
    int w, h;
	typename I::DataType value[pixelSize];
    if(ifs_in.is_open())
    {
        ifs_in >> w >> h;
        image.initItk(w, h);
        image.for_all_pixels([&] (typename I::PixelType &pix)
        {
			for(unsigned i=0; i<pixelSize; ++i)
				ifs_in >> value[i];
			memcpy(&pix, value, sizeof(typename I::PixelType));
        });
    }
    ifs_in.close();
	if(!ifs_in)
		return false;
	return true;
}

template <class I, class Compare>
void Histogram<I, Compare>::saveHistogram(const std::string& out, int nb_classes_per_dimension) const
{
    std::ofstream ofs_out(out);

    ofs_out << size() << std::endl;
    ofs_out << binsNumber() << std::endl;
    ofs_out << (nb_classes_per_dimension != 0 ? nb_classes_per_dimension : binsNumber()) << std::endl;

    for(const auto& bin : *this)
    {
        ofs_out << bin.first << " " << bin.second << std::endl;
    }
    ofs_out.close();

    return;
}

template <class I, class Compare>
void Histogram<I, Compare>::saveFullHistogram(const std::string& out) const
{
    std::ofstream ofs_out(out);

    ofs_out << size() << std::endl;
    ofs_out << binsNumber() << std::endl;
    for(const auto& bin : *this)
    {
        for(int i=0; i<bin.second; ++i)
            ofs_out << bin.first << std::endl;
    }
    ofs_out.close();

    return;
}

template <class I, class Compare>
void Histogram<I, Compare>::compute(const I& image)
{
    m_histogram.clear();

    m_size=image.width()*image.height();

    auto makeHisto = [this] (const PixelType& pix)
    {
        std::pair<iterator, bool> pos=m_histogram.insert(std::make_pair(pix, 0));
        ++(*pos.first).second;
    };

    image.for_all_pixels(makeHisto);

    return;
}

template <class I, class Compare>
void Histogram<I, Compare>::addPixel(const PixelType &bin, int occ)
{
    std::pair<iterator, bool> pos=m_histogram.insert(std::make_pair(bin, 0));
    (*pos.first).second+=occ;

    m_size += occ;
    return;
}

template <class I, class Compare>
void Histogram<I, Compare>::removePixel(const PixelType &bin, int occ)
{
    iterator pos=m_histogram.find(std::make_pair(bin, 0));
    if(pos!=m_histogram.end())
    {
        if(occ>=(*pos.first).second)
        {
            m_size -= (*pos.first).second;
            //also, remove the bin
            m_histogram.erase(pos);

        }
        else
        {
            m_size -= occ;
            (*pos.first).second -= occ;
        }
    }
    return;
}

template <class I, class Compare>
void Histogram<I, Compare>::setBin(const PixelType &bin, int occ)
{
    iterator pos=m_histogram.find(std::make_pair(bin, 0));
    if(occ==0)
    {
        //try to delete this bin
        if(pos!=m_histogram.end())
        {
            int previousOcc=(*pos).second;
            m_histogram.erase(pos);
            m_size-=previousOcc;
        }
    }
    else
    {
        int previousOcc= pos!=m_histogram.end() ? (*pos).second : 0;
        std::pair<iterator, bool> pos=m_histogram.insert(std::make_pair(bin, 0));
        (*pos.first).second=occ;

        m_size = m_size - previousOcc + occ;
    }
    return;
}

template <class I, class Compare>
void Histogram<I, Compare>::removeBin(const PixelType &bin)
{
    iterator pos=m_histogram.find(std::make_pair(bin, 0));
    if(pos!=m_histogram.end())
    {
        m_size -= (*pos.first).second;
        m_histogram.erase(pos);
    }
    return;
}

template <class I, class Compare>
void Histogram<I, Compare>::clear()
{
    m_histogram.clear();
    m_size=0;

    return;
}

template <class I, class Compare>
template<typename F_bothExist, typename F_firstExists, typename F_secondExists>
typename Histogram<I, Compare>::real Histogram<I, Compare>::compareB2BWith(const Histogram<I, Compare>& other,
                                        F_bothExist &lmbd_matchH1H2,
                                        F_firstExists &lmbd_matchH1,
                                        F_secondExists &lmbd_matchH2) const
{
    real d=0;                   //< distance

    if(&other==this)            //in case somebody has fun comparing H1 to H1
        return d;

    Compare lmbd_compare;             //< a class to compare the order of two bins
    PixelType p0, p1;           //< current bins
    int a0, a1;                 //< number of occurences for this bin
    real fa0, fa1;              //< represents a0 and a1 in frequency

    auto its = make_tuple(begin(), other.begin()); //< a double iterator over the histograms
    while(std::get<0>(its)!=end() && std::get<1>(its)!=other.end())
    {
        p0=(*std::get<0>(its)).first;
        a0=(*std::get<0>(its)).second;

        p1=(*std::get<1>(its)).first;
        a1=(*std::get<1>(its)).second;

        fa0=(real)a0/size();
        fa1=(real)a1/size();

        //compare p0 and p1
        if(lmbd_compare(p0, p1)) //no occurence of p0 in H1
        {
            d+=lmbd_matchH1(fa0);
            ++std::get<0>(its);
        }
        else if(lmbd_compare(p1, p0)) //no occurence of p1 in H0
        {
            d+=lmbd_matchH2(fa1);
            ++std::get<1>(its);
        }
        else //p0=p1 and they are both in H0 and H1
        {
            d+=lmbd_matchH1H2(fa0, fa1);
            ++std::get<0>(its);
            ++std::get<1>(its);
        }
    }

    while(std::get<0>(its)!=end()) //finish iteration over first histogram if necessary
    {
        a0=(*std::get<0>(its)).second;
        fa0=(real)a0/size();
        d+=lmbd_matchH1(fa0);
        ++std::get<0>(its);
    }
    while(std::get<1>(its)!=other.end()) //finish iteration over second histogram if necessary
    {
        a1=(*std::get<1>(its)).second;
        fa1=(real)a1/size();
        d+=lmbd_matchH2(fa1);
        ++std::get<1>(its);
    }

    return d;
}

template <class I, class Compare>
typename Histogram<I, Compare>::real Histogram<I, Compare>::compareIntersectionWith(const Histogram<I, Compare>& other) const
{
    int sizeOther = other.size();

    auto lmbd_intersection_match = [sizeOther](real fa0, real fa1) -> real
    {
        return std::min(fa0, fa1)/sizeOther;
    };

    auto lmbd_intersection_noMatch = [](real) -> real
    {
        return 0;
    };

    return 1.0 - compareB2BWith(other, lmbd_intersection_match, lmbd_intersection_noMatch, lmbd_intersection_noMatch);
}

template <class I, class Compare>
template <int L_norm>
typename Histogram<I, Compare>::real Histogram<I, Compare>::compareMinkowskiWith(const Histogram<I, Compare>& other) const
{
    auto lmbd_minkowski_match = [](real fa0, real fa1) -> real
    {
        return std::pow((fa0-fa1), (real)L_norm);
    };

    auto lmbd_minkowski_noMatch = [](real f) -> real
    {
        return std::pow(f, (real)L_norm);
    };

    return std::pow(compareB2BWith(other, lmbd_minkowski_match, lmbd_minkowski_noMatch, lmbd_minkowski_noMatch), 1.0/L_norm) / size();
}

template <class I, class Compare>
typename Histogram<I, Compare>::real Histogram<I, Compare>::compareChi2With(const Histogram<I, Compare>& other) const
{
    auto lmbd_chi2_match = [](real fa0, real fa1) -> real
    {
        real m_i = (fa0+fa1)/2;
        return (fa0-m_i)*(fa1-m_i)/m_i;
    };

    auto lmbd_chi2_matchH1 = [](real f) -> real
    {
        real m_i = f/2;
        return (f-m_i)*(f-m_i)/m_i;
    };

    auto lmbd_chi2_matchH2 = [](real f) -> real
    {
        real m_i = f/2;
        return m_i;
    };

    return compareB2BWith(other, lmbd_chi2_match, lmbd_chi2_matchH1, lmbd_chi2_matchH2);
}

template <class I, class Compare>
template <typename NUMBER, typename std::enable_if<std::is_arithmetic<NUMBER>::value>::type*>
typename Histogram<I, Compare>::real Histogram<I, Compare>::fitsUniformKS(NUMBER inf, NUMBER sup) const
{
    real d=0;         //< our distance
    int i=0;            //< at what sample we're at

    auto scaleBetween0and1 = [inf, sup] (NUMBER x) -> real //< used to normalize between 0 and 1
    {
        return ((real)x - inf) / (sup - inf);
    };

    for(const auto& it : *this)
    {
        //check https://www.cs.indiana.edu/~kapadia/project2/node14.html for the formula
        //note how our histogram is already sorted using Compare(p0, p1)

        real x=scaleBetween0and1(it.first);
        i+=1;
        d = std::max( d, std::abs(x - (real)(i-1)/size()) );  //< if we look at the first element of the bin, we want to compare it with D-
        i+=it.second-1;
        d = std::max( d, std::abs(x - ((real)i/size())) );    //< if we look at the last element of the bin, we want to compare it with D+
    }

    return d;
}

//
//Class to compare RGB pixels as less<itkRGB> does not produce the intended effect
//

template <typename T>
class CompareRGBPixels_lexicographic
{
public:
    CompareRGBPixels_lexicographic() {}
    bool operator()(const typename ImageCommon<ImageRGBBase<T>, false>::PixelType& object, const typename ImageCommon<ImageRGBBase<T>, false>::PixelType& other) const
    {
        int i;
        for(i=0; object[i] == other[i] && i<2; ++i);
        return object[i] < other[i];
    }
};

//
//Specialized histogram for RGB images
//

template <typename T>
class HistogramRGBBase : public Histogram<ImageCommon<ImageRGBBase<T>, false>, CompareRGBPixels_lexicographic<T>>
{
public:

    typedef typename Histogram<ImageCommon<ImageRGBBase<T>, false>, CompareRGBPixels_lexicographic<T>>::real real;
    typedef typename Histogram<ImageCommon<ImageRGBBase<T>, false>, CompareRGBPixels_lexicographic<T>>::PixelType PixelType;

    typedef CompareRGBPixels_lexicographic<T> Compare;

    HistogramRGBBase();
    HistogramRGBBase(const ImageCommon<ImageRGBBase<T>, false>& image);

    /**
     * @brief quantize turns a histogram into a clusterized integer histogram.
     * @param inf minimum of the pixelType (any value bellow is considered at minimum)
     * @param sup maximum of the pixelType (any value above is considered at maximum)
     * @param nb_classes_per_dimension the number of clusters (classes) per dimension.
     * 0 means it will choose automatically.
     * @return The clusterized histogram.
     */
    HistogramRGBBase<int> quantize(PixelType inf, PixelType sup, int nb_classes_per_dimension=0) const;

    /**
     * @brief updateStatistics updates the statistics mean and covariance.
     * Generally called automatically, one should still call it after changing manually the histogram data.
     */
    void updateStatistics();

    /**
     * @brief compute gets a histogram from an image.
     * Called automatically after a construction, one should still call it after using the default constructor.
     * @param image
     */
    void compute(const ImageCommon<ImageRGBBase<T>, false>& image);


    const real           mean(int i) const {return m_mean[i];}
    const real&    covariance(int i, int j) const {if(i==j) return m_covariance[i]; else return m_covariance[2 + i + j];}

    void clear();

    PixelType meanPixelType() const;

private:

    real      m_mean[3];
    real      m_covariance[6];
};

template <typename T>
HistogramRGBBase<T>::HistogramRGBBase() :
    Histogram<ImageCommon<ImageRGBBase<T>, false>, CompareRGBPixels_lexicographic<T>>(), m_mean(), m_covariance()
{
}

template <typename T>
HistogramRGBBase<T>::HistogramRGBBase(const ImageCommon<ImageRGBBase<T>, false>& image) :
    Histogram<ImageCommon<ImageRGBBase<T>, false>, CompareRGBPixels_lexicographic<T>>(image), m_mean(), m_covariance()
{
    compute(image);
}

template <typename T>
void HistogramRGBBase<T>::compute(const ImageCommon<ImageRGBBase<T>, false>& image)
{
    Histogram<ImageCommon<ImageRGBBase<T>, false>, CompareRGBPixels_lexicographic<T>>::compute(image);

    updateStatistics();

    return;
}

template <typename T>
void HistogramRGBBase<T>::updateStatistics()
{
    for(int i=0; i<3; ++i)
        m_mean[i] = 0;
    for(int i=0; i<6; ++i)
        m_covariance[i] = 0;

    for(const auto& bin : *this)
    {
        for(int i=0; i<3; ++i)
        {
            m_mean[i] += (real)bin.first[i]*bin.second;
        }
    }
    for(int i=0; i<3; ++i)
        m_mean[i] /= this->size();

    //here we have a mean vector so we can compute the covariance
    for(const auto& bin : *this)
    {
        for(int i=0; i<3; ++i)
        {
            real deviation=bin.first[i]-m_mean[i];
            for(int j=0; j<3; ++j)
            {
                //cov. matrix is set like this : /[a x x]     0 1 2 3 4 5
                                                //[A b x] -> [a b c A B C] (compact representation)
                                                //[B C c]
                if(i==j) //lower case (diagonal)
                    m_covariance[i] += bin.second * deviation * (bin.first[j]-m_mean[j]);
                else if(i>j) //upper case (others)
                    m_covariance[2 + i + j] += bin.second * deviation * (bin.first[j]-m_mean[j]);
                //ignore case x (repetition)
            }
        }

    }

    //m_covariance /= (this->size()-1);
    // v

    for(int i=0; i<6; ++i)
        m_covariance[i] /= (this->size()-1);

    return;
}

template <typename T>
HistogramRGBBase<int> HistogramRGBBase<T>::quantize(PixelType inf, PixelType sup, int nb_classes_per_dimension) const
{
    HistogramRGBBase<int> histogramQuantized;

    //first, compute the number of classes per dimension if it is not explicitely provided (nb_classes_per_dimension=0) or wrong (<0)

    if(nb_classes_per_dimension<=0)
    {
        nb_classes_per_dimension=(int)std::ceil(std::log2(this->size())+1.0); //Sturges formula
    }

    //fill histogram
    for(const auto& bin : *this)
    {
        int intervalInt3[3];
        for(int i=0; i<3; ++i)
        {
            intervalInt3[i] = std::max(0, std::min(nb_classes_per_dimension-1, (int) ( (real)(bin.first[i] - inf[i]) * nb_classes_per_dimension / (sup[i] - inf[i]) )));
        }
        HistogramRGBBase<int>::PixelType currentInterval( intervalInt3 );
        histogramQuantized.addPixel(currentInterval, bin.second);
    }

    //update mean
    histogramQuantized.updateStatistics();
    return histogramQuantized;
}

template <typename T>
void HistogramRGBBase<T>::clear()
{
    Histogram<ImageCommon<ImageRGBBase<T>, false>, CompareRGBPixels_lexicographic<T>>::clear();
    for(int i=0; i<3; ++i)
    {
        m_mean[i] = 0;
    }
    //there should be no need to clear covariance, but remember it won't be cleared

    return;
}

template <typename T>
typename HistogramRGBBase<T>::PixelType HistogramRGBBase<T>::meanPixelType() const
{
    PixelType mean;
    for(int i=0; i<3; ++i)
    {
        mean[i]=(T)m_mean[i];

    }

    return mean;
}

//
//Specialized histogram for Gray images
//

template <typename T>
class HistogramGrayBase : public Histogram<ImageCommon<ImageGrayBase<T>, false>, std::less<T>>
{
public:
    HistogramGrayBase();
    HistogramGrayBase(const ImageCommon<ImageGrayBase<T>, false>& image);

    typedef typename Histogram<ImageCommon<ImageGrayBase<T>, false>, std::less<T>>::real real;
    typedef typename Histogram<ImageCommon<ImageGrayBase<T>, false>, std::less<T>>::PixelType PixelType;

    void compute(const ImageCommon<ImageGrayBase<T>, false>& image);

    //avoid scalar character types to be displayed as such
    void saveHistogram(const std::string& out, int nb_classes_per_dimension=0) const;
    void saveFullHistogram(const std::string& out) const;

    real            min() {return m_min;}
    real            max() {return m_max;}
    real            mean() {return m_mean;}
    real            variance() {return m_variance;}

    void updateStatistics();

    void            clear();

    /**
     * @brief meanPixelType
     * @return the mean except it's of the type of the PixelType. (why again?)
     */
    PixelType       meanPixelType();

    /**
     * @brief quantize turns a histogram into a clusterized integer histogram.
     * @param inf minimum of the pixelType (any value bellow is considered at minimum)
     * @param sup maximum of the pixelType (any value above is considered at maximum)
     * @param nb_classes the number of clusters (classes). 0 means it will choose automatically.
     * @return The clusterized histogram.
     */
    HistogramGrayBase<int> quantize(PixelType inf=PixelType(0), PixelType sup=PixelType(1), unsigned int nb_classes=0) const;

    //test of fit

    real            fitsNormalChi2() const;

private:

    real            m_min;
    real            m_max;
    real            m_mean;
    real            m_variance;
};

template <typename T>
HistogramGrayBase<T>::HistogramGrayBase() :
    Histogram<ImageCommon<ImageGrayBase<T>, false>, std::less<T>>(),
    m_min(std::numeric_limits<int>::max()),
    m_max(std::numeric_limits<int>::min()),
    m_mean(0),
    m_variance(0)
{}

template <typename T>
HistogramGrayBase<T>::HistogramGrayBase(const ImageCommon<ImageGrayBase<T>, false>& image) :
    Histogram<ImageCommon<ImageGrayBase<T>, false>, std::less<T>>(image), m_mean(0), m_variance(0)
{
    compute(image);
}

template <typename T>
void HistogramGrayBase<T>::compute(const ImageCommon<ImageGrayBase<T>, false>& image)
{
    Histogram<ImageCommon<ImageGrayBase<T>, false>, std::less<T>>::compute(image);

    updateStatistics();

    return;
}

template <typename T>
void HistogramGrayBase<T>::saveHistogram(const std::string& out, int nb_classes_per_dimension) const
{
    std::ofstream ofs_out(out);

    ofs_out << this->size() << std::endl;
    ofs_out << this->binsNumber() << std::endl;
    ofs_out << (nb_classes_per_dimension != 0 ? nb_classes_per_dimension : this->binsNumber()) << std::endl;
    for(const auto& bin : *this)
    {
        ofs_out << +bin.first << " " << bin.second << std::endl;
    }
    ofs_out.close();

    return;
}

template <typename T>
void HistogramGrayBase<T>::saveFullHistogram(const std::string& out) const
{
    std::ofstream ofs_out(out);

    ofs_out << this->size() << std::endl;
    ofs_out << this->binsNumber() << std::endl;
    for(const auto& bin : *this)
    {
        for(int i=0; i<bin.second; ++i)
            ofs_out << +bin.first << std::endl; //hacked to print numbers no matter what type T is
    }
    ofs_out.close();

    return;
}

template <typename T>
void HistogramGrayBase<T>::updateStatistics()
{
    //reset
    m_min=std::numeric_limits<int>::max();
    m_max=std::numeric_limits<int>::min();
    m_mean=0;
    m_variance=0;
    //mean, min and max
    for(const auto& bin : *this)
    {
        m_min = std::min(m_min, real(bin.first));
        m_max = std::max(m_max, real(bin.first));
        m_mean+=(real)bin.first*bin.second / this->size();
    }

    for(const auto& bin : *this)
    {
        real deviation=bin.first-m_mean;
        m_variance += bin.second*deviation*deviation;
    }
    m_variance/=(this->size()-1);
    return;
}

template <typename T>
void HistogramGrayBase<T>::clear()
{
    Histogram<ImageCommon<ImageGrayBase<T>, false>, std::less<T>>::clear();

    m_mean=0;
    m_variance=0;
    m_min=std::numeric_limits<int>::max();
    m_max=std::numeric_limits<int>::min();

    return;
}

template <typename T>
typename HistogramGrayBase<T>::PixelType HistogramGrayBase<T>::meanPixelType()
{
    PixelType mean;
    mean=(PixelType)m_mean;
    return mean;
}

template <typename T>
HistogramGrayBase<int> HistogramGrayBase<T>::quantize(PixelType inf, PixelType sup, unsigned int nb_classes) const
{
    HistogramGrayBase<int> histogramQuantized;

    //first, compute the number of classes if it is not explicitely provided (nb_classes=0)

    if(nb_classes==0)
    {
        nb_classes=(unsigned int)std::log(this->size())+2.0;
    }

    //then, fill the histogram
    for(const auto& bin : *this)
    {
        int currentInterval = (int) ( (real)(bin.first - inf) * nb_classes / (sup - inf) );
        histogramQuantized.addPixel(currentInterval, bin.second);
    }
    histogramQuantized.updateStatistics();

    return histogramQuantized;
}

template <typename T>
typename HistogramGrayBase<T>::real HistogramGrayBase<T>::fitsNormalChi2() const
{
    real d=0;
    for(const auto& bin : *this)
    {
        //to change
        real chiNotSquared=((real)(bin.first*bin.second)/this->size() - m_mean)/m_variance;
        d+=chiNotSquared*chiNotSquared;
    }
    return d;
}


//typedefs

typedef HistogramRGBBase<unsigned char>     HistogramRGBu8;
typedef HistogramRGBBase<double>            HistogramRGBd;
//etc...

typedef HistogramGrayBase<unsigned char>    HistogramGrayu8;
typedef HistogramGrayBase<double>           HistogramGrayd;
//etc...

typedef HistogramGrayBase<double>           HistogramSpectrald;
//etc...

//Class to compare histograms

#endif
