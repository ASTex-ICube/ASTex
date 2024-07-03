//
// Created by grenier on 12/10/23.
//

#ifndef CCVT_TEST_CHARLINE_ENRICHED_SEGMENT_H
#define CCVT_TEST_CHARLINE_ENRICHED_SEGMENT_H


template <class Segment, class Vertex_handle>
class CEnrichedSegment
{
private:
    Segment       m_segment;
    Vertex_handle m_left;
    Vertex_handle m_right;

public:
    CEnrichedSegment() { }

    CEnrichedSegment(const Segment& segment,
                     Vertex_handle left,
                     Vertex_handle right)
    {
        m_segment = segment;
        m_left  = left;
        m_right = right;
    }

    CEnrichedSegment& operator = (const CEnrichedSegment& rhs)
    {
        m_segment = rhs.m_segment;
        m_left  = rhs.m_left;
        m_right = rhs.m_right;
        return *this;
    }

    const Segment& segment() const
    {
        return m_segment;
    }

    Vertex_handle left() const
    {
        return m_left;
    }

    Vertex_handle right() const
    {
        return m_right;
    }
};


#endif //CCVT_TEST_CHARLINE_ENRICHED_SEGMENT_H
