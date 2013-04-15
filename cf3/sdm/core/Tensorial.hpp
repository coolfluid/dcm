// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef cf3_sdm_core_Tensorial_hpp
#define cf3_sdm_core_Tensorial_hpp

#include "cf3/common/BoostArray.hpp"
#include "cf3/common/StringConversion.hpp"

#include "cf3/math/MatrixTypes.hpp"
#include "cf3/math/Consts.hpp"

#include "cf3/sdm/core/ShapeFunction.hpp"
#include "cf3/sdm/core/LibCore.hpp"

#include <boost/assign/list_of.hpp>
using namespace boost::assign;

namespace cf3 {
namespace sdm {
namespace core {

////////////////////////////////////////////////////////////////////////////////

struct sdm_core_API Lagrange
{
  static Real coeff(const Real& ksi, const RealVector& pts, const Uint idx)
  {
    // declarations for efficiency
    Uint k;
    const Uint nb_pts(pts.size());

    Real prod=1;
    for(k=0; k<nb_pts;++k)
    {
      if (k!=idx)
        prod*= (ksi-pts[k])/(pts[idx]-pts[k]);
    }
    return prod;
  }

  static Real deriv_coeff(const Real& ksi, const RealVector& pts, const Uint idx)
  {
    // declarations for efficiency
    Uint t;
    Uint k;
    Real term;
    const Uint nb_pts(pts.size());

     Real deriv_coeff = 0.;
     for (t = 0; t < nb_pts; ++t)
     {
       if (t != idx)
       {
         term = 1./(pts[idx]-pts[t]);
         for (k = 0; k < nb_pts; ++k)
         {
           if (k != idx && k != t)
           {
             term *= (ksi-pts[k])/(pts[idx]-pts[k]);
           }
         }
         deriv_coeff += term;
       }
     }
     return deriv_coeff;
  }
};

////////////////////////////////////////////////////////////////////////////////

template< typename DISTRIBUTION_1D, int P >
class sdm_core_API Point : public sdm::core::ShapeFunction
{
public: // typedefs

  typedef boost::shared_ptr<Point>       Ptr;
  typedef boost::shared_ptr<Point const> ConstPtr;

public: // functions

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW  ///< storing fixed-sized Eigen structures

  static std::string type_name() { return "Point<"+common::to_str(P)+">"; }

  Point(const std::string& name = type_name())
    : sdm::core::ShapeFunction(name),
      m_flx_pt_dirs(1,KSI),
      m_order(P),
      m_zero(0.),
      m_one(1.)
  {
    m_sol_pts.resize(1,3);
    m_sol_pts << 0, 0, 0;
    m_flx_pts.resize(1,3);
    m_flx_pts << 0, 0, 0;

    m_face_flx_pts.resize(0);
    m_flx_pt_sign.resize(1,1.);
  }

  virtual ~Point() {}

  virtual void compute_value(const RealVector& local_coordinate, RealRowVector& value) const
  {
    value[0] = 1.;
  }
  virtual void compute_gradient(const RealVector& local_coordinate, RealMatrix& gradient) const
  {
    gradient.setZero();
  }
  virtual void compute_flux_value(const Uint orientation, const RealVector& local_coordinate, RealRowVector& value) const
  {
    value[0] = 1.;
  }
  virtual void compute_flux_derivative(const Uint orientation, const RealVector& local_coordinate, RealVector& derivative) const
  {
    derivative.setZero();
  }
 
  virtual mesh::GeoShape::Type shape() const { return mesh::GeoShape::POINT; }
  virtual Uint dimensionality() const { return DIM_3D; }
  virtual Uint nb_faces() const { return 0; }
  virtual Uint order() const { return m_order; }
  virtual Uint nb_sol_pts() const { return 1; }
  virtual Uint nb_flx_pts() const { return 1; }
  virtual const RealMatrix& sol_pts() const { return m_sol_pts; }
  virtual const RealMatrix& flx_pts() const { return m_flx_pts; }
  virtual Uint flx_pt_dir(const Uint flx_pt) const { cf3_assert(flx_pt<nb_flx_pts()); return m_flx_pt_dirs[flx_pt]; }
  virtual const RealMatrix& face_normals() const { return m_face_normals; }
  virtual const std::vector<Uint>& interior_flx_pts() const { return m_interior_flx_pts; }
  virtual const std::vector<Uint>& face_flx_pts(const Uint face_idx, const Uint orientation, const Uint rotation) const { cf3_assert(face_idx<nb_faces()); return m_face_flx_pts[face_idx]; }
  virtual const std::vector<Uint>& inverted_face_flx_pts(const Uint face_idx, const Uint rotation = 0) const { cf3_assert(face_idx<nb_faces()); return m_face_flx_pts[face_idx]; }
  virtual const Real& flx_pt_sign(const Uint flx_pt) const { return m_flx_pt_sign[flx_pt]; }

private: // data

  Real m_zero;
  Real m_one;

  Uint                                m_order;            ///< Order of the solution shape function
  RealMatrix                          m_flx_pts;          ///< Flux point coordinates
  RealMatrix                          m_sol_pts;          ///< Solution point coordinates
  std::vector<Uint>                   m_flx_pt_dirs;      ///< Per flux point, the directions this flux point contributes to
  RealMatrix                          m_face_normals;     ///< Rows are normals to faces according to FaceNumbering
  std::vector<Uint>                   m_interior_flx_pts; ///< Flux points that lie inside the cell, not on the faces
  std::vector<std::vector<Uint> >     m_face_flx_pts;     ///< Flux points that on the cell faces
  std::vector<Real>                   m_flx_pt_sign;      ///< Sign to be multiplied with computed flux in flx_pt in direction dir

};

////////////////////////////////////////////////////////////////////////////////

template< typename DISTRIBUTION_1D, int P >
class sdm_core_API Line : public sdm::core::ShapeFunction
{
public: // typedefs

  typedef boost::shared_ptr<Line>       Ptr;
  typedef boost::shared_ptr<Line const> ConstPtr;

private: // typedefs

  enum FaceNumbering {KSI_NEG=0, KSI_POS=1};

public: // functions

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW  ///< storing fixed-sized Eigen structures

  static std::string type_name() { return "Line<"+common::to_str(P)+">"; }

  Line(const std::string& name = type_name())
    : sdm::core::ShapeFunction(name),
      m_dist_1d( P ),
      m_order(m_dist_1d.nb_sol_pts-1)
  {
    m_sol_pts.resize(m_dist_1d.nb_sol_pts,1);
    for (Uint s=0; s<m_dist_1d.nb_sol_pts; ++s)
      m_sol_pts(s,KSI)=m_dist_1d.sol_pts[s];

    m_flx_pts.resize(m_dist_1d.nb_flx_pts,1);
    for (Uint f=0; f<m_dist_1d.nb_flx_pts; ++f)
      m_flx_pts(f,KSI)=m_dist_1d.flx_pts[f];

    m_face_flx_pts.resize(2);
    m_flx_pt_sign.resize(m_dist_1d.nb_flx_pts,1.);

    for (Uint f=0; f<m_dist_1d.nb_flx_pts; ++f)
    {
      if (f==0)
      {
        m_face_flx_pts[KSI_NEG].push_back(f);
        m_flx_pt_sign[f]=-1.;
      }
      else if(f==m_dist_1d.nb_flx_pts-1)
      {
        m_face_flx_pts[KSI_POS].push_back(f);
        m_flx_pt_sign[f]=+1.;
      }
      else
      {
        m_interior_flx_pts.push_back(f);
        m_flx_pt_sign[f]=+1.;
      }
    }

    m_face_normals.resize(nb_faces(),DIM_1D); m_face_normals.setZero();
    m_face_normals(KSI_NEG,KSI)=-1;
    m_face_normals(KSI_POS,KSI)=+1;

  }

  virtual ~Line() {}

  virtual void compute_value(const RealVector& local_coordinate, RealRowVector& value) const
  {
    cf3_assert(value.size()==m_dist_1d.nb_sol_pts);
    for (Uint s=0; s<m_dist_1d.nb_sol_pts; ++s) {
      value[s] = Lagrange::coeff(local_coordinate[KSI],m_dist_1d.sol_pts,s);
    }
  }
  virtual void compute_gradient(const RealVector& local_coordinate, RealMatrix& gradient) const
  {
    cf3_assert(gradient.rows()==DIM_1D);
    cf3_assert(gradient.cols()==m_dist_1d.nb_sol_pts);
    for (Uint s=0; s<m_dist_1d.nb_sol_pts; ++s) {
      gradient(KSI,s) = Lagrange::deriv_coeff(local_coordinate[KSI],m_dist_1d.sol_pts,s);
    }
  }
  virtual void compute_flux_value(const Uint orientation, const RealVector& local_coordinate, RealRowVector& value) const
  {
    cf3_assert(value.size()==nb_flx_pts());
    cf3_assert(orientation==KSI);
    value.setZero();
    for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_flx_pts; ++f_ksi)
      value[f_ksi] = Lagrange::coeff(local_coordinate[KSI],m_dist_1d.flx_pts,f_ksi);
  }
  virtual void compute_flux_derivative(const Uint orientation, const RealVector& local_coordinate, RealVector& derivative) const
  {
    cf3_assert(derivative.size()==nb_flx_pts());
    cf3_assert(orientation==KSI);
    derivative.setZero();
    for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_flx_pts; ++f_ksi)
      derivative[f_ksi] = Lagrange::deriv_coeff(local_coordinate[KSI],m_dist_1d.flx_pts,f_ksi);
  }

  virtual mesh::GeoShape::Type shape() const { return mesh::GeoShape::LINE; }
  virtual Uint dimensionality() const { return DIM_1D; }
  virtual Uint nb_faces() const { return 2; }
  virtual Uint order() const { return m_order; }
  virtual Uint nb_sol_pts() const { return m_dist_1d.nb_sol_pts; }
  virtual Uint nb_flx_pts() const { return m_dist_1d.nb_flx_pts; }
  virtual const RealMatrix& sol_pts() const { return m_sol_pts; }
  virtual const RealMatrix& flx_pts() const { return m_flx_pts; }
  virtual Uint flx_pt_dir(const Uint flx_pt) const { cf3_assert(flx_pt<nb_flx_pts()); return KSI; }
  virtual const RealMatrix& face_normals() const { return m_face_normals; }
  virtual const std::vector<Uint>& interior_flx_pts() const { return m_interior_flx_pts; }
  virtual const std::vector<Uint>& face_flx_pts(const Uint face_idx, const Uint orientation, const Uint rotation) const { cf3_assert(face_idx<nb_faces()); return m_face_flx_pts[face_idx]; }
  virtual const std::vector<Uint>& inverted_face_flx_pts(const Uint face_idx, const Uint rotation = 0) const { cf3_assert(face_idx<nb_faces()); return m_face_flx_pts[face_idx]; }
  virtual const Real& flx_pt_sign(const Uint flx_pt) const { return m_flx_pt_sign[flx_pt]; }

private: // data

  DISTRIBUTION_1D                     m_dist_1d;          ///< holds 1D interpolation matrices
  Uint                                m_order;            ///< Order of the solution shape function
  RealMatrix                          m_flx_pts;          ///< Flux point coordinates
  RealMatrix                          m_sol_pts;          ///< Solution point coordinates
  RealMatrix                          m_face_normals;     ///< Rows are normals to faces according to FaceNumbering
  std::vector<Uint>                   m_interior_flx_pts; ///< Flux points that lie inside the cell, not on the faces
  std::vector< std::vector<Uint> >    m_face_flx_pts;     ///< Flux points that on the cell faces
  std::vector<Real>                   m_flx_pt_sign;                          ///< Sign to be multiplied with computed flux in flx_pt in direction dir

};

////////////////////////////////////////////////////////////////////////////////

template< typename DISTRIBUTION_1D, int P >
class sdm_core_API Quad : public sdm::core::ShapeFunction
{
public: // typedefs

  typedef boost::shared_ptr<Quad>       Ptr;
  typedef boost::shared_ptr<Quad const> ConstPtr;

private: // typedefs

  enum FaceNumbering {ETA_NEG=0, KSI_POS=1, ETA_POS=2, KSI_NEG=3};

  enum { NDIM=2 };
  enum { NFACES=4 };

public: // functions

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW  ///< storing fixed-sized Eigen structures

  static std::string type_name() { return "Quad<"+common::to_str(P)+">"; }


  Quad(const std::string& name = type_name())
    : sdm::core::ShapeFunction(name),
      m_dist_1d(P)
  {
    m_order = m_dist_1d.nb_sol_pts-1;
    m_nb_sol_pts = m_dist_1d.nb_sol_pts*m_dist_1d.nb_sol_pts;
    m_nb_flx_pts = m_dist_1d.nb_sol_pts*m_dist_1d.nb_flx_pts*DIM_2D;
    m_sol_pts.resize(m_nb_sol_pts,DIM_2D);
    m_flx_pts.resize(m_nb_flx_pts,DIM_2D);
    m_flx_pt_dirs.resize(m_nb_flx_pts);
    m_flx_pt_sign.resize(m_nb_flx_pts,1.);

    for (Uint s_ksi=0; s_ksi<m_dist_1d.nb_sol_pts; ++s_ksi)
    {
      for (Uint s_eta=0; s_eta<m_dist_1d.nb_sol_pts; ++s_eta)
      {
        const Uint s = s_eta*m_dist_1d.nb_sol_pts+s_ksi;

        m_sol_pts(s,KSI) = m_dist_1d.sol_pts[s_ksi];
        m_sol_pts(s,ETA) = m_dist_1d.sol_pts[s_eta];
      }
    }

    for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_flx_pts; ++f_ksi)
    {
      for (Uint f_eta=0; f_eta<m_dist_1d.nb_sol_pts; ++f_eta)
      {
        const Uint f = f_eta*m_dist_1d.nb_flx_pts+f_ksi;
        m_flx_pts(f,KSI) = m_dist_1d.flx_pts[f_ksi];
        m_flx_pts(f,ETA) = m_dist_1d.sol_pts[f_eta];
        m_flx_pt_dirs[f]=KSI;

        // f_eta is ignored as 1) the location may not be on faces; 2) it doesn't count as a face-point in the locally-1D line
        if (f_ksi==0)
        {
          m_flx_pt_sign[f]= -1.;
        }
        else if(f_ksi==m_dist_1d.nb_flx_pts-1)
        {
          m_flx_pt_sign[f]= +1.;
        }
        else
        {
          m_interior_flx_pts.push_back(f);
          m_flx_pt_sign[f]= +1.;
        }
      }
    }
    for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_sol_pts; ++f_ksi)
    {
      for (Uint f_eta=0; f_eta<m_dist_1d.nb_flx_pts; ++f_eta)
      {
        const Uint f = f_ksi*m_dist_1d.nb_flx_pts+f_eta + m_dist_1d.nb_sol_pts*m_dist_1d.nb_flx_pts;
        m_flx_pts(f,KSI) = m_dist_1d.sol_pts[f_ksi];
        m_flx_pts(f,ETA) = m_dist_1d.flx_pts[f_eta];
        m_flx_pt_dirs[f]=ETA;

        for (Uint s_eta=0; s_eta<m_dist_1d.nb_sol_pts; ++s_eta)
        {
          const Uint s = s_eta*m_dist_1d.nb_sol_pts + f_ksi;
        }

        // f_ksi is ignored as 1) the location may not be on faces; 2) it doesn't count as a face-point in the locally-1D line
        if (f_eta==0)
        {
          m_flx_pt_sign[f]=-1.;
        }
        else if(f_eta==m_dist_1d.nb_flx_pts-1)
        {
          m_flx_pt_sign[f]=+1.;
        }
        else
        {
          m_interior_flx_pts.push_back(f);
          m_flx_pt_sign[f]=+1.;
        }

      }
    }
    m_face_normals.resize(4,DIM_2D); m_face_normals.setZero();
    m_face_normals(KSI_NEG,KSI)=-1;
    m_face_normals(KSI_POS,KSI)=+1;
    m_face_normals(ETA_NEG,ETA)=-1;
    m_face_normals(ETA_POS,ETA)=+1;


    {
      Uint F=m_dist_1d.nb_flx_pts;
      Uint S=m_dist_1d.nb_sol_pts;

      Uint N_ksi=0;
      Uint N_eta=F*S;
      Uint Sm1=S-1;
      Uint Fm1=F-1;

      m_face_flx_pts.resize(NFACES, std::vector< std::vector<Uint> > (2, std::vector<Uint>(S)));

      for (Uint i=0; i<S; ++i)
      {
          m_face_flx_pts[KSI_NEG][0][i] = N_ksi      +  F*(Sm1-i);

          m_face_flx_pts[KSI_POS][0][i] = N_ksi+Fm1  +  F*i;

          m_face_flx_pts[ETA_NEG][0][i] = N_eta      +  F*i;

          m_face_flx_pts[ETA_POS][0][i] = N_eta+Fm1  +  F*(Sm1-i);

          m_face_flx_pts[KSI_NEG][1][i] = N_ksi      +  F*i;

          m_face_flx_pts[KSI_POS][1][i] = N_ksi+Fm1  +  F*(Sm1-i);

          m_face_flx_pts[ETA_NEG][1][i] = N_eta      +  F*(Sm1-i);

          m_face_flx_pts[ETA_POS][1][i] = N_eta+Fm1  +  F*i;

      }
    }
  }


  virtual ~Quad() {}
  virtual void compute_value(const RealVector& local_coordinate, RealRowVector& value) const
  {
    cf3_assert(value.size()==nb_sol_pts());
    for (Uint s_ksi=0; s_ksi<m_dist_1d.nb_sol_pts; ++s_ksi) {
      for (Uint s_eta=0; s_eta<m_dist_1d.nb_sol_pts; ++s_eta) {
        Uint s = s_eta*m_dist_1d.nb_sol_pts+s_ksi;
        value[s] = Lagrange::coeff(local_coordinate[KSI],m_dist_1d.sol_pts,s_ksi) * Lagrange::coeff(local_coordinate[ETA],m_dist_1d.sol_pts,s_eta);
      }
    }
  }
  virtual void compute_gradient(const RealVector& local_coordinate, RealMatrix& gradient) const
  {
    cf3_assert(gradient.rows()==DIM_2D);
    cf3_assert(gradient.cols()==nb_sol_pts());
    for (Uint s_ksi=0; s_ksi<m_dist_1d.nb_sol_pts; ++s_ksi) {
      for (Uint s_eta=0; s_eta<m_dist_1d.nb_sol_pts; ++s_eta) {
        Uint s = s_eta*m_dist_1d.nb_sol_pts+s_ksi;
        gradient(KSI,s) = Lagrange::deriv_coeff(local_coordinate[KSI],m_dist_1d.sol_pts,s_ksi) * Lagrange::coeff(local_coordinate[ETA],m_dist_1d.sol_pts,s_eta);
        gradient(ETA,s) = Lagrange::coeff(local_coordinate[KSI],m_dist_1d.sol_pts,s_ksi) * Lagrange::deriv_coeff(local_coordinate[ETA],m_dist_1d.sol_pts,s_eta);
      }
    }
  }
  virtual void compute_flux_value(const Uint orientation, const RealVector& local_coordinate, RealRowVector& value) const
  {
    cf3_assert(value.size()==nb_flx_pts());
    value.setZero();
    switch (orientation)
    {
    case KSI:
      for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_flx_pts; ++f_ksi) {
        for (Uint f_eta=0; f_eta<m_dist_1d.nb_sol_pts; ++f_eta) {
          const Uint f = f_eta*m_dist_1d.nb_flx_pts+f_ksi;
          value[f] = Lagrange::coeff(local_coordinate[KSI],m_dist_1d.flx_pts,f_ksi) * Lagrange::coeff(local_coordinate[ETA],m_dist_1d.sol_pts,f_eta);
        }
      }
      break;
    case ETA:
      for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_sol_pts; ++f_ksi) {
        for (Uint f_eta=0; f_eta<m_dist_1d.nb_flx_pts; ++f_eta)
        {
          const Uint f = f_ksi*m_dist_1d.nb_flx_pts+f_eta + m_dist_1d.nb_sol_pts*m_dist_1d.nb_flx_pts;
          value[f] = Lagrange::coeff(local_coordinate[KSI],m_dist_1d.sol_pts,f_ksi) * Lagrange::coeff(local_coordinate[ETA],m_dist_1d.flx_pts,f_eta);
        }
      }
      break;
    }
  }
  virtual void compute_flux_derivative(const Uint orientation, const RealVector& local_coordinate, RealVector& derivative) const
  {
    cf3_assert(derivative.size()==nb_flx_pts());
    derivative.setZero();
    switch (orientation)
    {
    case KSI:
      for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_flx_pts; ++f_ksi) {
        for (Uint f_eta=0; f_eta<m_dist_1d.nb_sol_pts; ++f_eta) {
          const Uint f = f_eta*m_dist_1d.nb_flx_pts+f_ksi;
          derivative[f] = Lagrange::deriv_coeff(local_coordinate[KSI],m_dist_1d.flx_pts,f_ksi) * Lagrange::coeff(local_coordinate[ETA],m_dist_1d.sol_pts,f_eta);
        }
      }
      break;
    case ETA:
      for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_sol_pts; ++f_ksi) {
        for (Uint f_eta=0; f_eta<m_dist_1d.nb_flx_pts; ++f_eta)
        {
          const Uint f = f_ksi*m_dist_1d.nb_flx_pts+f_eta + m_dist_1d.nb_sol_pts*m_dist_1d.nb_flx_pts;
          derivative[f] = Lagrange::coeff(local_coordinate[KSI],m_dist_1d.sol_pts,f_ksi) * Lagrange::deriv_coeff(local_coordinate[ETA],m_dist_1d.flx_pts,f_eta);
        }
      }
      break;
    }
  }

  virtual mesh::GeoShape::Type shape() const { return mesh::GeoShape::QUAD; }
  virtual Uint dimensionality() const { return DIM_2D; }
  virtual Uint nb_faces() const { return m_face_flx_pts.size(); }
  virtual Uint order() const { return m_order; }
  virtual Uint nb_sol_pts() const { return m_nb_sol_pts; }
  virtual Uint nb_flx_pts() const { return m_nb_flx_pts; }
  virtual const RealMatrix& sol_pts() const { return m_sol_pts; }
  virtual const RealMatrix& flx_pts() const { return m_flx_pts; }
  virtual Uint flx_pt_dir(const Uint flx_pt) const { cf3_assert(flx_pt<nb_flx_pts()); return m_flx_pt_dirs[flx_pt]; }
  virtual const RealMatrix& face_normals() const { return m_face_normals; }
  virtual const std::vector<Uint>& interior_flx_pts() const { return m_interior_flx_pts; }
  virtual const std::vector<Uint>& face_flx_pts(const Uint face_idx, const Uint orientation, const Uint rotation) const { cf3_assert(face_idx<nb_faces()); return m_face_flx_pts[face_idx][rotation]; }
  virtual const Real& flx_pt_sign(const Uint flx_pt) const { return m_flx_pt_sign[flx_pt]; }

private: // data

  DISTRIBUTION_1D                     m_dist_1d;                       ///< holds 1D interpolation matrices
  Uint                                m_order;                         ///< Order of the solution shape function
  Uint                                m_nb_sol_pts;                    ///< Number of solution points
  Uint                                m_nb_flx_pts;                    ///< Number of flux points
  RealMatrix                          m_flx_pts;                       ///< Flux point coordinates
  RealMatrix                          m_sol_pts;                       ///< Solution point coordinates
  std::vector<Uint>                   m_flx_pt_dirs;                   ///< Per flux point, the directions this flux point contributes to
  RealMatrix                          m_face_normals;                  ///< Rows are normals to faces according to FaceNumbering
  std::vector<Uint>                   m_interior_flx_pts;              ///< Flux points that lie inside the cell, not on the faces
  std::vector< std::vector< std::vector<Uint> > > m_face_flx_pts;      ///< Flux points that on the cell faces
  std::vector<Real>                   m_flx_pt_sign;                   ///< Sign to be multiplied with computed flux in flx_pt in direction dir
};


////////////////////////////////////////////////////////////////////////////////

template< typename DISTRIBUTION_1D, int P >
class sdm_core_API Hexa : public sdm::core::ShapeFunction
{
public: // typedefs

  typedef boost::shared_ptr<Hexa>       Ptr;
  typedef boost::shared_ptr<Hexa const> ConstPtr;

private: // typedefs

  enum FaceNumbering {ZTA_NEG=0, ZTA_POS=1, ETA_NEG=2, KSI_POS=3, ETA_POS=4, KSI_NEG=5};
  enum { NDIM=3 };
  enum { NFACES=6 };
  enum { NROTATIONS=4 };

public: // functions

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW  ///< storing fixed-sized Eigen structures

  static std::string type_name() { return "Hexa<"+common::to_str(P)+">"; }

  Hexa(const std::string& name = type_name())
    : sdm::core::ShapeFunction(name),
      m_dist_1d(P)
  {
    m_order = m_dist_1d.nb_sol_pts-1;
    m_nb_sol_pts = m_dist_1d.nb_sol_pts*m_dist_1d.nb_sol_pts*m_dist_1d.nb_sol_pts;
    m_nb_flx_pts = m_dist_1d.nb_sol_pts*m_dist_1d.nb_sol_pts*m_dist_1d.nb_flx_pts*DIM_3D;
    m_sol_pts.resize(m_nb_sol_pts,DIM_3D);
    m_flx_pts.resize(m_nb_flx_pts,DIM_3D);
    m_flx_pt_dirs.resize(m_nb_flx_pts);
    m_flx_pt_sign.resize(m_nb_flx_pts,1.);

    // Define solution points
    for (Uint s_ksi=0; s_ksi<m_dist_1d.nb_sol_pts; ++s_ksi)
    {
      for (Uint s_eta=0; s_eta<m_dist_1d.nb_sol_pts; ++s_eta)
      {
        for (Uint s_zta=0; s_zta<m_dist_1d.nb_sol_pts; ++s_zta)
        {
          const Uint s = (s_zta*m_dist_1d.nb_sol_pts + s_eta)*m_dist_1d.nb_sol_pts+s_ksi;

          m_sol_pts(s,KSI) = m_dist_1d.sol_pts[s_ksi];
          m_sol_pts(s,ETA) = m_dist_1d.sol_pts[s_eta];
          m_sol_pts(s,ZTA) = m_dist_1d.sol_pts[s_zta];
        }
      }
    }

    // define KSI direction flux points
    for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_flx_pts; ++f_ksi)
    {
      for (Uint f_eta=0; f_eta<m_dist_1d.nb_sol_pts; ++f_eta)
      {
        for (Uint f_zta=0; f_zta<m_dist_1d.nb_sol_pts; ++f_zta)
        {
          const Uint f = (f_zta*m_dist_1d.nb_sol_pts + f_eta)*m_dist_1d.nb_flx_pts+f_ksi;
          m_flx_pts(f,KSI) = m_dist_1d.flx_pts[f_ksi];
          m_flx_pts(f,ETA) = m_dist_1d.sol_pts[f_eta];
          m_flx_pts(f,ZTA) = m_dist_1d.sol_pts[f_zta];
          m_flx_pt_dirs[f]=KSI;

          // f_eta and f_zta are ignored as 1) the location may not be on faces; 2) it doesn't count as a face-point in the locally-1D line
          if (f_ksi==0)
          {
            m_flx_pt_sign[f]= -1.;
          }
          else if(f_ksi==m_dist_1d.nb_flx_pts-1)
          {
            m_flx_pt_sign[f]= +1.;
          }
          else
          {
            m_interior_flx_pts.push_back(f);
            m_flx_pt_sign[f]= +1.;
          }
        }
      }
    }
    // define ETA direction flux points
    for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_sol_pts; ++f_ksi)
    {
      for (Uint f_eta=0; f_eta<m_dist_1d.nb_flx_pts; ++f_eta)
      {
        for (Uint f_zta=0; f_zta<m_dist_1d.nb_sol_pts; ++f_zta)
        {
          const Uint f = (f_zta*m_dist_1d.nb_sol_pts + f_ksi)*m_dist_1d.nb_flx_pts+f_eta + m_dist_1d.nb_sol_pts*m_dist_1d.nb_sol_pts*m_dist_1d.nb_flx_pts;
          m_flx_pts(f,KSI) = m_dist_1d.sol_pts[f_ksi];
          m_flx_pts(f,ETA) = m_dist_1d.flx_pts[f_eta];
          m_flx_pts(f,ZTA) = m_dist_1d.sol_pts[f_zta];

          m_flx_pt_dirs[f]=ETA;

          // f_ksi and f_zta are ignored as 1) the location may not be on faces; 2) it doesn't count as a face-point in the locally-1D line
          if (f_eta==0)
          {
            m_flx_pt_sign[f]=-1.;
          }
          else if(f_eta==m_dist_1d.nb_flx_pts-1)
          {
            m_flx_pt_sign[f]=+1.;
          }
          else
          {
            m_interior_flx_pts.push_back(f);
            m_flx_pt_sign[f]=+1.;
          }

        }
      }
    }
    // define ZTA direction flux points
    for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_sol_pts; ++f_ksi)
    {
      for (Uint f_eta=0; f_eta<m_dist_1d.nb_sol_pts; ++f_eta)
      {
        for (Uint f_zta=0; f_zta<m_dist_1d.nb_flx_pts; ++f_zta)
        {
          const Uint f = (f_eta*m_dist_1d.nb_sol_pts + f_ksi)*m_dist_1d.nb_flx_pts+f_zta + 2*m_dist_1d.nb_sol_pts*m_dist_1d.nb_sol_pts*m_dist_1d.nb_flx_pts;
          m_flx_pts(f,KSI) = m_dist_1d.sol_pts[f_ksi];
          m_flx_pts(f,ETA) = m_dist_1d.sol_pts[f_eta];
          m_flx_pts(f,ZTA) = m_dist_1d.flx_pts[f_zta];

          m_flx_pt_dirs[f]=ZTA;

          // f_ksi and f_eta are ignored as 1) the location may not be on faces; 2) it doesn't count as a face-point in the locally-1D line
          if (f_zta==0)
          {
            m_flx_pt_sign[f]=-1.;
          }
          else if(f_zta==m_dist_1d.nb_flx_pts-1)
          {
            m_flx_pt_sign[f]=+1.;
          }
          else
          {
            m_interior_flx_pts.push_back(f);
            m_flx_pt_sign[f]=+1.;
          }

        }
      }
    }

    m_face_normals.resize(6,DIM_3D);
    // set all components to zero
    m_face_normals.setZero();
    // change 1 component for each face in the right direction
    m_face_normals(KSI_NEG,KSI)=-1;
    m_face_normals(KSI_POS,KSI)=+1;
    m_face_normals(ETA_NEG,ETA)=-1;
    m_face_normals(ETA_POS,ETA)=+1;
    m_face_normals(ZTA_NEG,ZTA)=-1;
    m_face_normals(ZTA_POS,ZTA)=+1;

    {
      Uint MATCHED=0;
      Uint INVERTED=1;
      Uint F=m_dist_1d.nb_flx_pts;
      Uint S=m_dist_1d.nb_sol_pts;

      Uint N_ksi=0;
      Uint N_eta=F*S*S;
      Uint N_zta=N_eta+F*S*S;
      Uint Sm1=S-1;
      Uint Fm1=F-1;

      m_face_flx_pts.resize(NFACES, std::vector< std::vector< std::vector<Uint> > >(2, std::vector< std::vector<Uint> >(NROTATIONS, std::vector<Uint>(S*S))));

      for (Uint i=0; i<S; ++i)
      {
        for (Uint j=0; j<S; ++j)
        {
          // The indices  0, 1, 2, 3  are the orientations (i.e. rotations of the cell compared to the face)
          m_face_flx_pts[KSI_NEG][MATCHED][0][i+S*j] = N_ksi      +  (F*S)*i        +  F*j;
          m_face_flx_pts[KSI_NEG][MATCHED][1][i+S*j] = N_ksi      +  F*i            +  (F*S)*(Sm1-j);
          m_face_flx_pts[KSI_NEG][MATCHED][2][i+S*j] = N_ksi      +  (F*S)*(Sm1-i)  +  F*(Sm1-j);
          m_face_flx_pts[KSI_NEG][MATCHED][3][i+S*j] = N_ksi      +  F*(Sm1-i)      +  (F*S)*j;

          m_face_flx_pts[KSI_POS][MATCHED][0][i+S*j] = N_ksi+Fm1  +  F*i            +  (F*S)*j;
          m_face_flx_pts[KSI_POS][MATCHED][1][i+S*j] = N_ksi+Fm1  +  (F*S)*i        +  F*(Sm1-j);
          m_face_flx_pts[KSI_POS][MATCHED][2][i+S*j] = N_ksi+Fm1  +  F*(Sm1-i)      +  (F*S)*(Sm1-j);
          m_face_flx_pts[KSI_POS][MATCHED][3][i+S*j] = N_ksi+Fm1  +  (F*S)*(Sm1-i)  +  F*j;

          m_face_flx_pts[ETA_NEG][MATCHED][0][i+S*j] = N_eta      +  F*i            +  (F*S)*j;
          m_face_flx_pts[ETA_NEG][MATCHED][1][i+S*j] = N_eta      +  (F*S)*i        +  F*(Sm1-j);
          m_face_flx_pts[ETA_NEG][MATCHED][2][i+S*j] = N_eta      +  F*(Sm1-i)      +  (F*S)*(Sm1-j);
          m_face_flx_pts[ETA_NEG][MATCHED][3][i+S*j] = N_eta      +  (F*S)*(Sm1-i)  +  F*j;

          m_face_flx_pts[ETA_POS][MATCHED][0][i+S*j] = N_eta+Fm1  +  (F*S)*i        +  F*j;
          m_face_flx_pts[ETA_POS][MATCHED][1][i+S*j] = N_eta+Fm1  +  F*i            +  (F*S)*(Sm1-j);
          m_face_flx_pts[ETA_POS][MATCHED][2][i+S*j] = N_eta+Fm1  +  (F*S)*(Sm1-i)  +  F*(Sm1-j);
          m_face_flx_pts[ETA_POS][MATCHED][3][i+S*j] = N_eta+Fm1  +  F*(Sm1-i)      +  (F*S)*j;

          m_face_flx_pts[ZTA_NEG][MATCHED][0][i+S*j] = N_zta      +  (F*S)*i        +  F*j;
          m_face_flx_pts[ZTA_NEG][MATCHED][1][i+S*j] = N_zta      +  F*i            +  (F*S)*(Sm1-j);
          m_face_flx_pts[ZTA_NEG][MATCHED][2][i+S*j] = N_zta      +  (F*S)*(Sm1-i)  +  F*(Sm1-j);
          m_face_flx_pts[ZTA_NEG][MATCHED][3][i+S*j] = N_zta      +  F*(Sm1-i)      +  (F*S)*j;

          m_face_flx_pts[ZTA_POS][MATCHED][0][i+S*j] = N_zta+Fm1  +  F*i            +  (F*S)*j;
          m_face_flx_pts[ZTA_POS][MATCHED][1][i+S*j] = N_zta+Fm1  +  (F*S)*i        +  F*(Sm1-j);
          m_face_flx_pts[ZTA_POS][MATCHED][2][i+S*j] = N_zta+Fm1  +  F*(Sm1-i)      +  (F*S)*(Sm1-j);
          m_face_flx_pts[ZTA_POS][MATCHED][3][i+S*j] = N_zta+Fm1  +  (F*S)*(Sm1-i)  +  F*j;

          m_face_flx_pts[KSI_NEG][INVERTED][0][i+S*j] = N_ksi      +  F*i            +  (F*S)*j;
          m_face_flx_pts[KSI_NEG][INVERTED][1][i+S*j] = N_ksi      +  (F*S)*(Sm1-i)  +  F*j;
          m_face_flx_pts[KSI_NEG][INVERTED][2][i+S*j] = N_ksi      +  F*(Sm1-i)      +  (F*S)*(Sm1-j);
          m_face_flx_pts[KSI_NEG][INVERTED][3][i+S*j] = N_ksi      +  (F*S)*i        +  F*(Sm1-j);

          m_face_flx_pts[KSI_POS][INVERTED][0][i+S*j] = N_ksi+Fm1  +  (F*S)*i        +  F*j;
          m_face_flx_pts[KSI_POS][INVERTED][1][i+S*j] = N_ksi+Fm1  +  F*(Sm1-i)      +  (F*S)*j;
          m_face_flx_pts[KSI_POS][INVERTED][2][i+S*j] = N_ksi+Fm1  +  (F*S)*(Sm1-i)  +  F*(Sm1-j);
          m_face_flx_pts[KSI_POS][INVERTED][3][i+S*j] = N_ksi+Fm1  +  F*i            +  (F*S)*(Sm1-j);

          m_face_flx_pts[ETA_NEG][INVERTED][0][i+S*j] = N_eta      +  (F*S)*i        +  F*j;
          m_face_flx_pts[ETA_NEG][INVERTED][1][i+S*j] = N_eta      +  F*(Sm1-i)      +  (F*S)*j;
          m_face_flx_pts[ETA_NEG][INVERTED][2][i+S*j] = N_eta      +  (F*S)*(Sm1-i)  +  F*(Sm1-j);
          m_face_flx_pts[ETA_NEG][INVERTED][3][i+S*j] = N_eta      +  F*i            +  (F*S)*(Sm1-j);

          m_face_flx_pts[ETA_POS][INVERTED][0][i+S*j] = N_eta+Fm1  +  F*i            +  (F*S)*j;
          m_face_flx_pts[ETA_POS][INVERTED][1][i+S*j] = N_eta+Fm1  +  (F*S)*(Sm1-i)  +  F*j;
          m_face_flx_pts[ETA_POS][INVERTED][2][i+S*j] = N_eta+Fm1  +  F*(Sm1-i)      +  (F*S)*(Sm1-j);
          m_face_flx_pts[ETA_POS][INVERTED][3][i+S*j] = N_eta+Fm1  +  (F*S)*i        +  F*(Sm1-j);

          m_face_flx_pts[ZTA_NEG][INVERTED][0][i+S*j] = N_zta      +  F*i            +  (F*S)*j;
          m_face_flx_pts[ZTA_NEG][INVERTED][1][i+S*j] = N_zta      +  (F*S)*(Sm1-i)  +  F*j;
          m_face_flx_pts[ZTA_NEG][INVERTED][2][i+S*j] = N_zta      +  F*(Sm1-i)      +  (F*S)*(Sm1-j);
          m_face_flx_pts[ZTA_NEG][INVERTED][3][i+S*j] = N_zta      +  (F*S)*i        +  F*(Sm1-j);

          m_face_flx_pts[ZTA_POS][INVERTED][0][i+S*j] = N_zta+Fm1  +  (F*S)*i        +  F*j;
          m_face_flx_pts[ZTA_POS][INVERTED][1][i+S*j] = N_zta+Fm1  +  F*(Sm1-i)      +  (F*S)*j;
          m_face_flx_pts[ZTA_POS][INVERTED][2][i+S*j] = N_zta+Fm1  +  (F*S)*(Sm1-i)  +  F*(Sm1-j);
          m_face_flx_pts[ZTA_POS][INVERTED][3][i+S*j] = N_zta+Fm1  +  F*i            +  (F*S)*(Sm1-j);

        }
      }
    }


  }


  virtual ~Hexa() {}
  virtual void compute_value(const RealVector& local_coordinate, RealRowVector& value) const
  {
    cf3_assert(value.size()==nb_sol_pts());
    for (Uint s_ksi=0; s_ksi<m_dist_1d.nb_sol_pts; ++s_ksi) {
      for (Uint s_eta=0; s_eta<m_dist_1d.nb_sol_pts; ++s_eta) {
        for (Uint s_zta=0; s_zta<m_dist_1d.nb_sol_pts; ++s_zta)
        {
          const Uint s = (s_zta*m_dist_1d.nb_sol_pts+s_eta)*m_dist_1d.nb_sol_pts+s_ksi;
          value[s] = Lagrange::coeff(local_coordinate[KSI],m_dist_1d.sol_pts,s_ksi)
                   * Lagrange::coeff(local_coordinate[ETA],m_dist_1d.sol_pts,s_eta)
                   * Lagrange::coeff(local_coordinate[ZTA],m_dist_1d.sol_pts,s_zta);
        }
      }
    }
  }
  virtual void compute_gradient(const RealVector& local_coordinate, RealMatrix& gradient) const
  {
    cf3_assert(gradient.rows()==DIM_3D);
    cf3_assert(gradient.cols()==nb_sol_pts());
    for (Uint s_ksi=0; s_ksi<m_dist_1d.nb_sol_pts; ++s_ksi) {
      for (Uint s_eta=0; s_eta<m_dist_1d.nb_sol_pts; ++s_eta) {
        for (Uint s_zta=0; s_zta<m_dist_1d.nb_sol_pts; ++s_zta)
        {
          const Uint s = (s_zta*m_dist_1d.nb_sol_pts+s_eta)*m_dist_1d.nb_sol_pts+s_ksi;
          gradient(KSI,s) = Lagrange::deriv_coeff(local_coordinate[KSI],m_dist_1d.sol_pts,s_ksi)
                          * Lagrange::coeff(local_coordinate[ETA],m_dist_1d.sol_pts,s_eta)
                          * Lagrange::coeff(local_coordinate[ZTA],m_dist_1d.sol_pts,s_zta);
          gradient(ETA,s) = Lagrange::coeff(local_coordinate[KSI],m_dist_1d.sol_pts,s_ksi)
                          * Lagrange::deriv_coeff(local_coordinate[ETA],m_dist_1d.sol_pts,s_eta)
                          * Lagrange::coeff(local_coordinate[ZTA],m_dist_1d.sol_pts,s_zta);
          gradient(ZTA,s) = Lagrange::coeff(local_coordinate[KSI],m_dist_1d.sol_pts,s_ksi)
                          * Lagrange::coeff(local_coordinate[ETA],m_dist_1d.sol_pts,s_eta)
                          * Lagrange::deriv_coeff(local_coordinate[ZTA],m_dist_1d.sol_pts,s_zta);
        }
      }
    }
  }
  virtual void compute_flux_value(const Uint orientation, const RealVector& local_coordinate, RealRowVector& value) const
  {
    cf3_assert(value.size()==nb_flx_pts());
    value.setZero();
    switch (orientation)
    {
    case KSI:
      for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_flx_pts; ++f_ksi) {
        for (Uint f_eta=0; f_eta<m_dist_1d.nb_sol_pts; ++f_eta) {
          for (Uint f_zta=0; f_zta<m_dist_1d.nb_sol_pts; ++f_zta) {
            const Uint f = (f_zta*m_dist_1d.nb_sol_pts+f_eta)*m_dist_1d.nb_flx_pts+f_ksi;
            value[f] = Lagrange::coeff(local_coordinate[KSI],m_dist_1d.flx_pts,f_ksi)
                     * Lagrange::coeff(local_coordinate[ETA],m_dist_1d.sol_pts,f_eta)
                     * Lagrange::coeff(local_coordinate[ZTA],m_dist_1d.sol_pts,f_zta);
          }
        }
      }
      break;
    case ETA:
      for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_sol_pts; ++f_ksi) {
        for (Uint f_eta=0; f_eta<m_dist_1d.nb_flx_pts; ++f_eta) {
          for (Uint f_zta=0; f_zta<m_dist_1d.nb_sol_pts; ++f_zta)
          {
            const Uint f = (f_zta*m_dist_1d.nb_sol_pts + f_ksi)*m_dist_1d.nb_flx_pts+f_eta + m_dist_1d.nb_sol_pts*m_dist_1d.nb_sol_pts*m_dist_1d.nb_flx_pts;
            value[f] = Lagrange::coeff(local_coordinate[KSI],m_dist_1d.sol_pts,f_ksi)
                     * Lagrange::coeff(local_coordinate[ETA],m_dist_1d.flx_pts,f_eta)
                     * Lagrange::coeff(local_coordinate[ZTA],m_dist_1d.sol_pts,f_zta);
          }
        }
      }
      break;
    case ZTA:
      for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_sol_pts; ++f_ksi) {
        for (Uint f_eta=0; f_eta<m_dist_1d.nb_sol_pts; ++f_eta) {
          for (Uint f_zta=0; f_zta<m_dist_1d.nb_flx_pts; ++f_zta)
          {
            const Uint f = (f_eta*m_dist_1d.nb_sol_pts + f_ksi)*m_dist_1d.nb_flx_pts+f_zta + 2*m_dist_1d.nb_sol_pts*m_dist_1d.nb_sol_pts*m_dist_1d.nb_flx_pts;
            value[f] = Lagrange::coeff(local_coordinate[KSI],m_dist_1d.sol_pts,f_ksi)
                     * Lagrange::coeff(local_coordinate[ETA],m_dist_1d.sol_pts,f_eta)
                     * Lagrange::coeff(local_coordinate[ZTA],m_dist_1d.flx_pts,f_zta);
          }
        }
      }
      break;
    }
  }
  virtual void compute_flux_derivative(const Uint orientation, const RealVector& local_coordinate, RealVector& derivative) const
  {
    cf3_assert(derivative.size()==nb_flx_pts());
    derivative.setZero();

    switch (orientation)
    {
    case KSI:
      for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_flx_pts; ++f_ksi) {
        for (Uint f_eta=0; f_eta<m_dist_1d.nb_sol_pts; ++f_eta) {
          for (Uint f_zta=0; f_zta<m_dist_1d.nb_sol_pts; ++f_zta) {
            const Uint f = (f_zta*m_dist_1d.nb_sol_pts+f_eta)*m_dist_1d.nb_flx_pts+f_ksi;
            derivative[f] = Lagrange::deriv_coeff(local_coordinate[KSI],m_dist_1d.flx_pts,f_ksi)
                          * Lagrange::coeff(local_coordinate[ETA],m_dist_1d.sol_pts,f_eta)
                          * Lagrange::coeff(local_coordinate[ZTA],m_dist_1d.sol_pts,f_zta);
          }
        }
      }
      break;
    case ETA:
      for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_sol_pts; ++f_ksi) {
        for (Uint f_eta=0; f_eta<m_dist_1d.nb_flx_pts; ++f_eta) {
          for (Uint f_zta=0; f_zta<m_dist_1d.nb_sol_pts; ++f_zta)
          {
            const Uint f = (f_zta*m_dist_1d.nb_sol_pts + f_ksi)*m_dist_1d.nb_flx_pts+f_eta + m_dist_1d.nb_sol_pts*m_dist_1d.nb_sol_pts*m_dist_1d.nb_flx_pts;
            derivative[f] = Lagrange::coeff(local_coordinate[KSI],m_dist_1d.sol_pts,f_ksi)
                          * Lagrange::deriv_coeff(local_coordinate[ETA],m_dist_1d.flx_pts,f_eta)
                          * Lagrange::coeff(local_coordinate[ZTA],m_dist_1d.sol_pts,f_zta);
          }
        }
      }
      break;
    case ZTA:
      for (Uint f_ksi=0; f_ksi<m_dist_1d.nb_sol_pts; ++f_ksi) {
        for (Uint f_eta=0; f_eta<m_dist_1d.nb_sol_pts; ++f_eta) {
          for (Uint f_zta=0; f_zta<m_dist_1d.nb_flx_pts; ++f_zta)
          {
            const Uint f = (f_eta*m_dist_1d.nb_sol_pts + f_ksi)*m_dist_1d.nb_flx_pts+f_zta + 2*m_dist_1d.nb_sol_pts*m_dist_1d.nb_sol_pts*m_dist_1d.nb_flx_pts;
            derivative[f] = Lagrange::coeff(local_coordinate[KSI],m_dist_1d.sol_pts,f_ksi)
                          * Lagrange::coeff(local_coordinate[ETA],m_dist_1d.sol_pts,f_eta)
                          * Lagrange::deriv_coeff(local_coordinate[ZTA],m_dist_1d.flx_pts,f_zta);
          }
        }
      }
      break;
    }

  }

  virtual mesh::GeoShape::Type shape() const { return mesh::GeoShape::HEXA; }
  virtual Uint dimensionality() const { return DIM_3D; }
  virtual Uint nb_faces() const { return NFACES; }
  virtual Uint order() const { return m_order; }
  virtual Uint nb_sol_pts() const { return m_nb_sol_pts; }
  virtual Uint nb_flx_pts() const { return m_nb_flx_pts; }
  virtual const RealMatrix& sol_pts() const { return m_sol_pts; }
  virtual const RealMatrix& flx_pts() const { return m_flx_pts; }
  virtual Uint flx_pt_dir(const Uint flx_pt) const { cf3_assert(flx_pt<nb_flx_pts()); return m_flx_pt_dirs[flx_pt]; }
  virtual const RealMatrix& face_normals() const { return m_face_normals; }
  virtual const std::vector<Uint>& interior_flx_pts() const { return m_interior_flx_pts; }
  virtual const std::vector<Uint>& face_flx_pts(const Uint face_idx, const Uint orientation, const Uint rotation) const { cf3_assert(face_idx<nb_faces()); cf3_assert(rotation<NROTATIONS); return m_face_flx_pts[face_idx][orientation][rotation]; }
  virtual const Real& flx_pt_sign(const Uint flx_pt) const { return m_flx_pt_sign[flx_pt]; }

private: // data

  DISTRIBUTION_1D                     m_dist_1d;                              ///< holds 1D interpolation matrices
  Uint                                m_order;                                ///< Order of the solution shape function
  Uint                                m_nb_sol_pts;                           ///< Number of solution points
  Uint                                m_nb_flx_pts;                           ///< Number of flux points
  RealMatrix                          m_flx_pts;                              ///< Flux point coordinates
  RealMatrix                          m_sol_pts;                              ///< Solution point coordinates
  std::vector< Uint >                 m_flx_pt_dirs;                          ///< Per flux point, the directions this flux point contributes to
  std::vector< std::vector< std::vector< std::vector<Uint> > > > m_face_flx_pts;      ///< Flux points that on the cell faces
  RealMatrix                          m_face_normals;                         ///< Rows are normals to faces according to FaceNumbering
  std::vector<Uint>                   m_interior_flx_pts;                     ///< Flux points that lie inside the cell, not on the faces
  std::vector<Real>                   m_flx_pt_sign;                          ///< Sign to be multiplied with computed flux in flx_pt in direction dir
};


////////////////////////////////////////////////////////////////////////////////

} // core
} // sdm
} // cf3

#endif // cf3_sdm_core_ShapeFunctionT_hpp
