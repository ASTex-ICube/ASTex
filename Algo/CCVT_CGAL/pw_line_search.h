//
// Created by grenier on 12/10/23.
//

#ifndef CCVT_TEST_CHARLINE_PW_LINE_SEARCH_H
#define CCVT_TEST_CHARLINE_PW_LINE_SEARCH_H

#include <iostream>

#include "line_search.h"


//------------//
// CWCWeights //
//------------//

template <class CCVT, class T>
class CLSWeights : public CLineSearch<T, T>
{
protected:
    CCVT* m_ccvt;
    unsigned m_nb;

public:
    CLSWeights(CCVT* ccvt,
               const unsigned max_iters,
               const double max_alpha)
            : CLineSearch<T, T>(max_iters, max_alpha)
    {
        m_ccvt = ccvt;
        m_nb = m_ccvt->count_visible_sites();
    }

    double compute_function() const
    {
        return m_ccvt->compute_wcvt_energy();
    }

    void compute_gradient(std::vector<T>& V) const
    {
        m_ccvt->compute_weight_gradient(V);
    }

    bool update_scene(const std::vector<T>& X)
    {
        if (m_ccvt->connectivity_fixed())
        {
            m_ccvt->clean_pixels();
            m_ccvt->assign_pixels();
            return true;
        }

        m_ccvt->update_weights(X, false);
        return m_ccvt->update_triangulation(true);
        //return has_same_vertices();
    }

    bool has_same_vertices() const
    {
        unsigned nb = m_ccvt->count_visible_sites();
        if (nb != m_nb) std::cout << "HiddenVertices: " << m_nb << " -> " << nb << std::endl;
        return (nb == m_nb);
    }
};

//-------------//
// CWCPosition //
//-------------//

template <class CCVT, class Position, class Velocity>
class CLSPositions : public CLineSearch<Position, Velocity>
{
protected:
    CCVT* m_ccvt;
    unsigned m_nb;

public:
    CLSPositions(CCVT* ccvt,
                 const unsigned max_iters,
                 const double max_alpha)
            : CLineSearch<Position, Velocity>(max_iters, max_alpha)
    {
        m_ccvt = ccvt;
        m_nb = m_ccvt->count_visible_sites();
    }

    double compute_function() const
    {
        return ( - m_ccvt->compute_wcvt_energy() );
    }

    void compute_gradient(std::vector<Velocity>& V) const
    {
        m_ccvt->compute_position_gradient(V, -1.0);
    }

    bool update_scene(const std::vector<Position>& X)
    {
        if (m_ccvt->connectivity_fixed())
        {
            m_ccvt->clean_pixels();
            m_ccvt->assign_pixels();
            return true;
        }

        m_ccvt->update_positions(X, true, false);
        return m_ccvt->update_triangulation(true);
        //return has_same_vertices();
    }

    bool has_same_vertices() const
    {
        unsigned nb = m_ccvt->count_visible_sites();
        if (nb != m_nb) std::cout << "HiddenVertices: " << m_nb << " -> " << nb << std::endl;
        return (nb == m_nb);
    }
};

#endif //CCVT_TEST_CHARLINE_PW_LINE_SEARCH_H
