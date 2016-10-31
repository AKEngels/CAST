#ifndef SCON_LINKEDCELLS_HEADER
#define SCON_LINKEDCELLS_HEADER

#include <vector>
#include <stddef.h>
#include <stdexcept>
#include <limits>
#include <utility>
#include "scon.h"
#include "scon_vect.h"
#include "scon_utility.h"
#include "configuration.h"
#pragma once


namespace scon
{

  template<class T3D>
  inline T3D clip_to_periodic_box(T3D const & value, T3D const & box)
  {
    return T3D(clip_to_circular_range(value.x(), -(box.x() / 2.), box.x() / 2.),
               clip_to_circular_range(value.y(), -(box.y() / 2.), box.y() / 2.),
               clip_to_circular_range(value.z(), -(box.z() / 2.), box.z() / 2.));
  }

  namespace linked
  {

    //template<class T, class V3D = c3<T>, class NV3D = std::vector<V3D>> class Cells;

    struct fragmentation { enum T { full = 1, half = 2 }; };

    // get box number of pos @param1
    template<class int_type = std::ptrdiff_t, class int3d_type = c3<int_type>>
    int_type box_number_from_offset(int3d_type const &i, int3d_type const & d, 
                                    bool const & periodics = false)
    {
      if (periodics)
      {
        int_type x = i.x() % d.x(), 
          y = i.y() % d.y(), 
          z = i.z() % d.z();
        return z*d.x()*d.y() + y*d.x() + x;
      }
      return i.z()*d.x()*d.y() + i.y()*d.x() + i.x();
    }

    template<class Box>
    class box_neighbour_iter
      : public std::iterator<std::bidirectional_iterator_tag, 
      const typename Box::int_type>
    {

      Box const * m_box;
      typename Box::int3d_type d;
      typename Box::int_type offset;
      bool m_end;

      void inc()
      {
        auto const & f = m_box->cells().fragments();
        if (d.x() < f) ++d.x();
        else
        {
          if (d.y() < f)
          {
            d.x() = -f;
            ++d.y();
          }
          else
          {
            if (d.z() < f)
            {
              d.x() = d.y() = -f;
              ++d.z();
            }
            else m_end = true;
          }
        }
      }

      bool dec()
      {
        auto const & f = m_box->cells().fragments();
        if (d.x() > -f) --d.x();
        else
        {
          if (d.y() > -f)
          {
            --d.y();
            d.x() = f;
          }
          else
          {
            if (d.z() > -f)
            {
              --d.z();
              d.y() = d.x() = f;
            }
            else return false;
          }
        }
        return true;
      }

      void inc_to_system()
      {
        while (!m_end && !in_system()) inc();
      }

      void dec_to_system()
      {
        bool deced = true;
        while (!in_system() && deced) {
          deced = dec();
        }
      }

      bool in_system()
      {
        auto q = d + m_box->offset();
        auto dq = (q - m_box->cells().dimensions()) + 1;
        if (q.x() < 0 || dq.x() > 0 || q.y() < 0 ||
          dq.y() > 0 || q.z() < 0 || dq.z() > 0)
        {
          return false;
        }
        return true;
      }

    public:

      box_neighbour_iter(Box const & box, bool const end = true)
        : m_box(&box), d(end ? m_box->cells().fragments() : -m_box->cells().fragments()),
        offset(-1), m_end(end)
      { 
        if (!m_box->cells().periodic())
        {
          if (m_end) dec_to_system();
          else inc_to_system();
        }
        offset = m_box->relative_box(d);
      }

      typedef box_neighbour_iter<Box>                      my_type;
      typedef std::bidirectional_iterator_tag              iterator_category;
      typedef std::iterator<iterator_category,
        const typename Box::int_type>                      iterator_type;
      typedef typename iterator_type::difference_type      difference_type;
      typedef typename iterator_type::value_type           value_type;
      typedef typename std::remove_const<value_type>::type non_const_value_type;
      typedef const non_const_value_type                   const_value_type;
      typedef typename iterator_type::pointer              pointer;
      typedef typename iterator_type::reference            reference;

      reference operator* () const { return offset; }
      pointer   operator->() const { return &offset; }

      bool operator== (my_type const &rhs) const 
      { 
        if (m_end == true && rhs.m_end == true) return true;
        if (m_end != rhs.m_end) return false;
        return (d == rhs.d && offset == rhs.offset); 
      }
      bool operator!= (my_type const &rhs) const { return !operator==(rhs); }
      bool operator> (my_type const &rhs) const
      {
        if (*this == rhs) return false;
        if (d.z() > rhs.d.z()) return true;
        if (d.z() == rhs.d.z() && d.y() > rhs.d.y()) return true;
        if (d.z() == rhs.d.z() && d.y() == rhs.d.y() && d.x() > rhs.d.x()) return true;
        return false;
      }
      bool operator< (my_type const &rhs) const
      {
        if (*this == rhs) return false;
        if (d.z() < rhs.d.z()) return true;
        if (d.z() == rhs.d.z() && d.y() < rhs.d.y()) return true;
        if (d.z() == rhs.d.z() && d.y() == rhs.d.y() && d.x() < rhs.d.x()) return true;
        return false;
      }
      bool operator<= (my_type const &rhs) const { return !operator>(rhs); }
      bool operator>= (my_type const &rhs) const { return !operator<(rhs); }

      typename Box::int3d_type const & delta() const { return d; }

      /*
      ++ increase priority => x > y > z
      x,y,z sequence
      0,0,0 -> 1,0,0 -> 2,0,0 -> 0,1,0 -> 1,1,0 -> 2,1,0 ...
      */
      my_type & operator++ ()
      {
        inc();
        if (!m_box->cells().periodic()) inc_to_system();
        offset = m_box->relative_box(d);
        return *this;
      }
      my_type & operator-- ()
      {
        auto const f = m_box->cells().fragments();
        if (m_end)
        {
          m_end = false;
        }
        else
        {
          if (!dec())
          {
            throw std::logic_error("Cannot iterate lower than -F(ragmentation), -F , -F in linked cells.");
          }
          if (!m_box->cells().periodic()) dec_to_system();
        }
        offset = m_box->relative_box(d);
        return *this;
      }
    };


    template<class Box>
    class box_element_iter
      : public std::iterator<std::forward_iterator_tag, 
      const typename Box::int_type>
    {
    public:
      typedef box_element_iter<Box>                        my_type;
      typedef std::forward_iterator_tag                    iterator_category;
      typedef std::iterator<iterator_category,
        const typename Box::int_type>                      iterator_type;
      typedef typename iterator_type::difference_type      difference_type;
      typedef typename iterator_type::value_type           value_type;
      typedef typename std::remove_const<value_type>::type non_const_value_type;
      typedef const non_const_value_type                   const_value_type;
      typedef typename iterator_type::pointer              pointer;
      typedef typename iterator_type::reference            reference;
    private:
      Box const * m_box;
      typename Box::int_type m_element;
    public:
      box_element_iter(Box const &box, typename Box::int_type element)
        : m_box(&box), m_element(element) {}
      bool operator== (box_element_iter const &rhs)
      {
        return (m_element == rhs.m_element && m_box == rhs.m_box);
      }
      bool operator!= (box_element_iter const &rhs)
      {
        return !operator==(rhs);
      }
      my_type& operator++ ()
      {
        if (m_element >= 0)
        {
          m_element = m_box->links()[static_cast<std::size_t>(m_element)];
        }
      }
      my_type& operator++ (int)
      {
        my_type ret(*this);
        ++(*this);
        return ret;
      }
      reference operator* () const { return m_element; }
      pointer   operator->() const { return &m_element; }
    };


    template<class Box>
    class box_adjacency_iter
      : public std::iterator<std::forward_iterator_tag, 
      const typename Box::int_type>
    {

    public:

      typedef box_adjacency_iter<Box>                      my_type;
      typedef std::forward_iterator_tag                    iterator_category;
      typedef std::iterator<iterator_category,
        const typename Box::int_type>                      iterator_type;
      typedef typename iterator_type::difference_type      difference_type;
      typedef typename iterator_type::value_type           value_type;
      typedef typename std::remove_const<value_type>::type non_const_value_type;
      typedef const non_const_value_type                   const_value_type;
      typedef typename iterator_type::pointer              pointer;
      typedef typename iterator_type::reference            reference;
      typedef typename Box::neighbour_iterator             neighbour_iterator;

    private:

      Box const * m_box;
      neighbour_iterator m_neighbour;
      neighbour_iterator m_nbe;
      typename Box::int_type m_element;
      bool m_end;

    public:

      box_adjacency_iter(Box const &box, bool const begin = false)
        : m_box(&box), m_neighbour(m_box->neighbours().begin()), 
        m_nbe(m_box->neighbours().end()),
        m_element(-1), m_end(true)
      {
        if (begin && m_neighbour != m_nbe)
        {
          auto cbox = *m_neighbour;
          //std::cout << "Adjacency begin iter for box " << cbox << "\n";
          if (cbox >= 0 && std::size_t(cbox) < m_box->cells().roots().size() &&
              m_box->cells().roots()[std::size_t(cbox)] >= 0)
          {
            m_element = m_box->cells().roots()[std::size_t(cbox)];
            //std::cout << "First valid element is " << m_element << "\n";
            m_end = false;
          }
          while (m_element < 0 && ++m_neighbour != m_nbe)
          {
            //std::cout << "Trying box " << *m_neighbour << "\n";
            if (*m_neighbour < 0) continue;
            std::size_t const sbox(static_cast<std::size_t>(*m_neighbour));
            // We're past the last box -> end reached
            if (sbox > m_box->cells().roots().size()) break;
            if (m_box->cells().roots()[sbox] >= 0)
            {
              // first valid element found
              m_element = m_box->cells().roots()[sbox];
              //std::cout << "First valid element is " << m_element << "\n";
              m_end = false;
              break;
            }
          }
        }
      }

      bool operator== (my_type const &rhs) const
      {
        if (m_end == true && rhs.m_end == true) return true;
        if (m_end != rhs.m_end) return false;
        return ((m_box == rhs.m_box) && m_neighbour == rhs.m_neighbour && m_element == rhs.m_element);
      }

      bool operator!= (my_type const &rhs) const 
      { 
        return !operator==(rhs);
      }

      reference operator* () const { return m_element; }
      pointer   operator->() const { return &m_element; }

      // increment and decrement
      my_type & operator++ (void)
      {
        if (m_end) return *this;
        m_end = true;
        // next atom in current box
        m_element = m_box->cells().links()[m_element];
        // end iterator
        //std::cout << "Increasing adja, new atom " << m_element << "\n";
        // if next_atom < 0 we have no more atoms in this box
        while (m_element < 0 && ++m_neighbour != m_nbe)
        {
          if (*m_neighbour < 0) continue;
          std::size_t const sbox((std::size_t)*m_neighbour);
          if (sbox >= m_box->cells().roots().size()) break;
          if (m_box->cells().roots()[sbox] >= 0)
          {
            m_element = m_box->cells().roots()[sbox];
            break;
          }
        }
        if (m_element >= 0 && std::size_t(m_element) < m_box->cells().size())
        {
          m_end = false;
        }
        //std::cout << " is end = " << m_end << "\n";
        return *this;
      }

      my_type operator++ (int)
      { // return copy then increment
        my_type ret(*this);
        ++(*this);
        return ret;
      }

    };

    template<class int_type = std::ptrdiff_t>
    inline c3<int_type> box_offset_by_id(int_type boxid, c3<int_type> const & dim)
    {
      auto const plane = dim.x()*dim.y();
      auto const z = boxid / plane;
      auto const rest = boxid % plane;
      auto const y = rest / dim.x();
      auto const x = rest % dim.x();
      return c3<int_type>{ x, y, z };
      //int_type const plane = dim.x()*dim.y();
      //int_type const modplane = boxid % plane;
      //int_type const modlane = modplane % dim.x();
      //return int3d_type(modlane, 
      //                  (modplane - modlane) / dim.x(), 
      //                  (boxid - modplane) / plane);
    }


    template<class Box>
    class box_proxy
    {
    public:

      typedef typename Box::int3d_type      int3d_type;
      typedef typename Box::int_type        int_type;
      typedef typename Box::vec_type        vec_type;
      typedef Box                           cells_type;
      typedef box_proxy<Box>                my_type;
      typedef box_neighbour_iter<my_type>   neighbour_iterator;
      typedef box_element_iter<my_type>     element_iterator;
      typedef box_adjacency_iter<my_type>   adjacent_iterator;

    private:

      Box const & m_lc;
      int3d_type const box_idx;
      int_type const m_index;
      box_proxy& operator= (box_proxy const &);

    public:

      box_proxy(Box const & links, std::size_t index)
        : m_lc(links), box_idx(box_offset_by_id<std::size_t>(index, scon::c3<std::size_t>{links.dimensions()})),
        m_index(index < std::size_t(std::numeric_limits<int_type>::max())?int_type(index):-1)
      { }

      int_type relative_box(int3d_type d) const
      {
        //std::cout << "Relative box of " << d;
        d += box_idx;
        //std::cout << " + " << box_idx << " = " << d << " in " << m_lc.dimensions();
        if (m_lc.periodic())
        {
          int_type x = d.x() % m_lc.dimensions().x(),
                   y = d.y() % m_lc.dimensions().y(), 
                   z = d.z() % m_lc.dimensions().z();
          if (x < 0) x += m_lc.dimensions().x();
          if (y < 0) y += m_lc.dimensions().y();
          if (z < 0) z += m_lc.dimensions().z();
          return z*m_lc.dimensions().x()*m_lc.dimensions().y() 
            + y*m_lc.dimensions().x() + x;
        }
        int_type ret;
        if (d.x() < 0 || d.y() < 0 || d.z() < 0 || 
          d.x() >= m_lc.dimensions().x() || 
          d.y() >= m_lc.dimensions().y() ||
          d.z() >= m_lc.dimensions().z())
        {
          ret = -1;
        }
        else
        {
          auto line = m_lc.dimensions().x();
          auto plane = line*m_lc.dimensions().y();
          ret = d.z()*plane + d.y()*line + d.x();
        }
        //std::cout << " is " << ret << "\n";
        return ret;
      }

      std::size_t id() const { return m_index; }

      cells_type const & cells() const
      {
        return m_lc;
      }

      vec_type position() const
      {
        return m_lc.box_position(m_index);
      }

      int3d_type offset() const
      {
        return m_lc.box_offset(m_index);
      }

      class neighbours_type
      {
        my_type const & m_box;
        neighbours_type& operator= (neighbours_type const &);
      public:
        typedef box_neighbour_iter<my_type> iterator;
        neighbours_type(my_type const &box) : m_box(box) { }
        iterator begin() const
        {
          return iterator(m_box, false);
        }
        iterator end() const
        {
          return iterator(m_box, true);
        }
      };

      neighbours_type neighbours() const
      {
        return neighbours_type(*this);
      }

      class adjacent_type
      {
        my_type const & m_box;
        adjacent_type& operator= (adjacent_type const &);
      public:
        typedef box_adjacency_iter<my_type> iterator;
        adjacent_type(my_type const &box) : m_box(box) { }
        iterator begin() const
        {
          return iterator(m_box, true);
        }
        iterator end() const
        {
          return iterator(m_box, false);
        }
      };

      adjacent_type adjacencies() const
      {
        return adjacent_type(*this);
      }

      int_type root() const
      {
        return m_lc.roots()[std::size_t(m_index)];
      }
      bool empty() const
      {
        return root() < 0;
      }
      element_iterator begin() const
      {
        return element_iterator(*this, root());
      }
      element_iterator end() const
      {
        return element_iterator(*this, -1);
      }
      std::vector<int_type> const & links() { return m_lc.links(); }
      bool operator== (my_type const &rhs) 
      { 
        return &m_lc == &rhs.m_lc && m_index == rhs.m_index; 
      }
      bool operator!= (my_type const &rhs) { return !operator==(rhs); }
    };



    template<class T, class V3D = c3<T>, 
	  class NV3D = scon::vector<V3D>>
    class Cells
    {
    public:

      typedef Cells<T, V3D, NV3D>      my_type;
      typedef V3D                      vec_type;
      typedef NV3D                     nvec_type;
      typedef std::ptrdiff_t           int_type;
      typedef scon::c3<std::ptrdiff_t> int3d_type;
      typedef box_proxy<my_type>       box_type;

    private:

      Cells& operator= (Cells const &);

      vec_type pb_clip(vec_type const &v) const
      {
        return clip_to_periodic_box<vec_type>(v, m_pb);
      }

      void setMaxMin(void)
      {
        using std::min;
        using std::max;
        m_min = vec_type(0);
        m_max = vec_type(0);
        if (!empty())
        {
          m_min = m_periodic ? pb_clip(positions.front()) : positions.front();
          m_max = m_min;
          for (auto const & position : positions)
          {
            vec_type const p = m_periodic ? 
              pb_clip(position) : position;
            m_min = min(p, m_min);
            m_max = max(p, m_max);
          }
          if (m_periodic)
          {
            vec_type const phb = m_pb 
              / T(2), mhb(-phb);
            m_min = min(mhb, m_min);
            m_max = max(phb, m_max);
          }
          // space extension
          m_max += extent;
          m_min -= extent;
          if (Config::get().general.verbosity > 29U)
          {
            std::cout << "LinkedCells::setMaxMin found max: ";
            std::cout << m_max << ", min: " << m_min << "." << endl;
          }
        }
      }

      void calcDimensions(void)
      {
        using std::floor;
        T const edgeInverse = 1.0 / edge;
        zero_diff = floor(m_min*edgeInverse);
        if (Config::get().general.verbosity > 29U)
        {
          std::cout << "LinkedCells::calcDimensions edgeinv: ";
          std::cout << edgeInverse << ", zerodiff: " << zero_diff << "." << endl;
        }
        auto dim_t = (floor(m_max*edgeInverse) - zero_diff) + 1.0;
        dim = int3d_type(std::ptrdiff_t(dim_t.x()), 
          std::ptrdiff_t(dim_t.y()), 
          std::ptrdiff_t(dim_t.z()));
        std::size_t const cells(static_cast<std::size_t>(dim.x()*dim.y()*dim.z()));
        if (Config::get().general.verbosity > 29U)
        {
          std::cout << "cells: " << cells << ", dim: " << dim << "." << endl;
        }
        m_roots.assign(cells, -1);
        Nbox.assign(cells, 0);
      }

      void link()
      {
        std::size_t const N(positions.size());
        for (std::size_t i(0u); i<N; ++i)
        {
          std::size_t const cell = box_index_from_point(positions[i]);
          if (Config::get().general.verbosity > 29U)
          {
            std::cout << "LinkedCells::link(): Box of element ";
            std::cout << i << " is " << cell << "." << endl;
          }
          cellofelement[i] = cell;
          m_links[i] = m_roots.at(cell);
          m_roots[cell] = i;
          ++Nbox[cell];
        }
      }

    public:

      vec_type                 m_max, m_min, m_pb, zero_diff;
      int3d_type               dim;
      T                        edge, extent;
      std::vector<int_type>    m_roots, m_links;
      std::vector<std::size_t> cellofelement, Nbox;
      int_type                 m_fragmentation;
      nvec_type const &        positions;
      bool                     m_periodic;

      Cells(nvec_type const & rep, T const edge_length = 10.0, bool const periodicity = false,
            vec_type const & periodic_boundaries = vec_type(0), T const ext = 0.0,
            fragmentation::T const fragments = fragmentation::full)
        : m_pb(periodic_boundaries), edge(edge_length), extent(ext), m_fragmentation(fragments > 0?fragments:1),
        positions(rep), m_periodic(periodicity)
      {
        edge = edge_length > 0.0 ? edge_length 
          / static_cast<double>(m_fragmentation) : 10.0;
        update();
      }

      // extend the space further than max and min
      void extend(T e_val) { extent = e_val < 0.0 ? -e_val : e_val; }

      // set links
      void update()
      {
        if (!empty())
        {
          if (Config::get().general.verbosity > 29U)
          {
            std::cout << "LinkedCells::update for ";
            std::cout << positions.size() << " elements in boxes of " << edge << " and " << fragments() << " fragments.\n\n";
          }
          cellofelement.assign(positions.size(), 0u);
          m_links.assign(positions.size(), -1);
          setMaxMin();
          calcDimensions();
          link();
        }
      }

      // initialize the grid
      void init(T const edge_length = 10.0, 
                linked::fragmentation::T const f = 
                  linked::fragmentation::full)
      {
        m_fragmentation = int_type(f > 0 ? f : 1);
        edge = edge_length > 0.0 ? edge_length 
          / T(m_fragmentation) : 10.0;
        update();
      }

      // get box indices of pos @param1
      int3d_type box_offset(vec_type p) const
      {
        using scon::min;
        using scon::max;
        if (m_periodic) p = pb_clip(p);
        p = min(m_max, max(m_min, p));
        return int3d_type{ floor(p / edge) - zero_diff };
      }

      int3d_type box_offset(std::ptrdiff_t const box_index) const
      {
        return box_offset_by_id<std::ptrdiff_t>(box_index, dim);
      }

      vec_type box_position(std::size_t const box_index) const
      {
        return (vec_type(box_offset(box_index))*edge + zero_diff*edge);
      }

      vec_type position(std::size_t const element_index) const
      {
        return m_periodic ? pb_clip(positions[element_index]) 
          : positions[element_index];
      }

      // get box number of pos @param1
      std::size_t box_index_from_point(vec_type const & p) const 
      { 
        return box_index_from_offset(box_offset(p)); 
      }

      // get box number of pos index @param1
      std::size_t box_index_from_element(std::size_t const p) const 
      { 
        return cellofelement[p]; 
      }

      // get box number of pos @param1
      std::size_t box_index_from_offset(int3d_type const &i) const
      {
        int_type const N(i.z()*dim.x()*dim.y() + i.y()*dim.x() + i.x());
        return N < 0 ? 0u : (std::size_t)N;
      }

      bool is_in_cells(vec_type const & p)
      {
        if (m_periodic) return true;
        if (p.x() < m_min.x() || p.x() > m_max.x()) return false;
        else if (p.y() < m_min.y() || p.y() > m_max.y()) return false;
        else if (p.z() < m_min.z() || p.z() > m_max.z()) return false;
        return true;
      }

      std::size_t size() const { return positions.size(); }
      bool empty() const { return positions.empty(); }

      nvec_type const & pos(void) { return positions; }

      // Box Proxys

      box_type box(std::size_t const box_index) const
      { 
        return box_type(*this, box_index);
      }
      box_type box_of_element(std::size_t const element_index) const
      { 
        return box_type(*this, cellofelement[element_index]);
      }
      box_type box_of_point(vec_type const & p) const
      { 
        return box_type(*this, box_index_from_point(p));
      }

      std::vector<std::ptrdiff_t> const & links(void) const { return m_links; }
      std::vector<std::ptrdiff_t> const & roots(void) const { return m_roots; }
      int3d_type const & dimensions() const { return dim; }

      int_type fragments() const 
      { 
        return m_periodic ? m_fragmentation +1 : m_fragmentation;
      }
      bool periodic() const { return m_periodic; }

      bool verify() const
      {
        std::size_t const N = positions.size();
        if (Nbox.size() != m_roots.size()) return false;
        if (N != cellofelement.size() || N != m_links.size()) return false;
        for (std::size_t i = 0; i < N; ++i)
        {
          vec_type const p = m_periodic ? pb_clip(positions[i]) : positions[i];
          vec_type const box_p = box_of_element(i).position();
          if (p.x() < box_p.x() || p.x() > (box_p.x() + edge)) return false;
          if (p.y() < box_p.y() || p.y() > (box_p.y() + edge)) return false;
          if (p.z() < box_p.z() || p.z() > (box_p.z() + edge)) return false;
        }
        return true;
      }

    };


  }

}

#endif
