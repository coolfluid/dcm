// Copyright (C) 2010-2013 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE "Test module for cf3::sdm"

#include <boost/test/unit_test.hpp>
#include <boost/assign/list_of.hpp>
#include "cf3/common/Log.hpp"
#include "cf3/common/Core.hpp"
#include "cf3/common/Environment.hpp"
#include "cf3/common/OptionList.hpp"
#include "cf3/common/OSystem.hpp"
#include "cf3/common/OSystemLayer.hpp"
#include "cf3/common/StringConversion.hpp"
#include "cf3/math/Consts.hpp"
#include "cf3/math/MatrixTypes.hpp"
#include "cf3/math/Defs.hpp"
#include "cf3/sdm/core/Tensorial.hpp"
#include "cf3/sdm/core/LegendreGaussEnd.hpp"

using namespace boost::assign;
using namespace cf3;
using namespace cf3::common;
using namespace cf3::mesh;
using namespace cf3::sdm::core;

//////////////////////////////////////////////////////////////////////////////

#define CF3_CHECK_EQUAL(x1,x2)\
do {\
  Real fraction=200*math::Consts::eps(); \
  if (x2 == 0) \
    BOOST_CHECK_SMALL(x1,fraction); \
  else if (x1 == 0) \
      BOOST_CHECK_SMALL(x2,fraction); \
  else \
    BOOST_CHECK_CLOSE_FRACTION(x1,x2,fraction);\
} while(false)

//////////////////////////////////////////////////////////////////////////////

struct sdm_MPITests_Fixture
{
  /// common setup for each test case
  sdm_MPITests_Fixture()
  {
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~sdm_MPITests_Fixture()
  {
  }
  /// possibly common functions used on the tests below


  /// common values accessed by all tests goes here
  int    m_argc;
  char** m_argv;

};

////////////////////////////////////////////////////////////////////////////////

BOOST_FIXTURE_TEST_SUITE( sdm_solver_TestSuite, sdm_MPITests_Fixture )

//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( init_mpi )
{
#ifdef test_is_mpi
  PE::Comm::instance().init(m_argc,m_argv);
#endif
  Core::instance().environment().options().set("log_level",(Uint)INFO);
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( test_P2_quad )
{
  enum FaceNumbering {ETA_NEG=0, KSI_POS=1, ETA_POS=2, KSI_NEG=3};
  enum CellOrientation {MATCHED=0, INVERTED=1};

  boost::shared_ptr< Quad<LegendreGaussEnd,2> > sf = allocate_component< Quad<LegendreGaussEnd,2> >("sf");

  // MATCHED, ALSO MEANS NOT ROTATED
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,0)[0], 8  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,0)[1], 4  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,0)[2], 0  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,0)[0], 3  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,0)[1], 7  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,0)[2], 11 );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,0)[0], 12 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,0)[1], 16 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,0)[2], 20 );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,0)[0], 23 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,0)[1], 19 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,0)[2], 15 );

  // INVERTED, ALSO MEANS ROTATED
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,1)[0], 0  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,1)[1], 4  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,1)[2], 8  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,1)[0], 11 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,1)[1], 7  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,1)[2], 3  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,1)[0], 20 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,1)[1], 16 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,1)[2], 12 );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,1)[0], 15 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,1)[1], 19 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,1)[2], 23 );

}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( test_P2_hexa )
{
  enum FaceNumbering {ZTA_NEG=0, ZTA_POS=1, ETA_NEG=2, KSI_POS=3, ETA_POS=4, KSI_NEG=5};
  enum CellOrientation {MATCHED=0, INVERTED=1};
  boost::shared_ptr< Hexa<LegendreGaussEnd,2> > sf = allocate_component< Hexa<LegendreGaussEnd,2> >("sf");

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,0)[0], 0   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,0)[1], 12  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,0)[2], 24  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,0)[3], 4   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,0)[4], 16  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,0)[5], 28  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,0)[6], 8   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,0)[7], 20  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,0)[8], 32  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,1)[0], 24  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,1)[1], 28  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,1)[2], 32  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,1)[3], 12  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,1)[4], 16  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,1)[5], 20  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,1)[6], 0   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,1)[7], 4   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,1)[8], 8   );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,2)[0], 32  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,2)[1], 20  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,2)[2], 8   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,2)[3], 28  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,2)[4], 16  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,2)[5], 4   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,2)[6], 24  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,2)[7], 12  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,2)[8], 0   );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,3)[0], 8   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,3)[1], 4   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,3)[2], 0   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,3)[3], 20  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,3)[4], 16  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,3)[5], 12  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,3)[6], 32  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,3)[7], 28  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,MATCHED,3)[8], 24  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,0)[0], 3   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,0)[1], 7   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,0)[2], 11  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,0)[3], 15  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,0)[4], 19  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,0)[5], 23  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,0)[6], 27  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,0)[7], 31  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,0)[8], 35  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,1)[0], 11  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,1)[1], 23  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,1)[2], 35  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,1)[3], 7   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,1)[4], 19  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,1)[5], 31  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,1)[6], 3   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,1)[7], 15  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,1)[8], 27  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,2)[0], 35  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,2)[1], 31  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,2)[2], 27  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,2)[3], 23  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,2)[4], 19  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,2)[5], 15  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,2)[6], 11  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,2)[7], 7   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,2)[8], 3   );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,3)[0], 27  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,3)[1], 15  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,3)[2], 3   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,3)[3], 31  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,3)[4], 19  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,3)[5], 7   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,3)[6], 35  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,3)[7], 23  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,MATCHED,3)[8], 11  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,0)[0], 36  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,0)[1], 40  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,0)[2], 44  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,0)[3], 48  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,0)[4], 52  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,0)[5], 56  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,0)[6], 60  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,0)[7], 64  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,0)[8], 68  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,1)[0], 44  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,1)[1], 56  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,1)[2], 68  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,1)[3], 40  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,1)[4], 52  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,1)[5], 64  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,1)[6], 36  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,1)[7], 48  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,1)[8], 60  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,2)[0], 68  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,2)[1], 64  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,2)[2], 60  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,2)[3], 56  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,2)[4], 52  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,2)[5], 48  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,2)[6], 44  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,2)[7], 40  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,2)[8], 36  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,3)[0], 60  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,3)[1], 48  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,3)[2], 36  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,3)[3], 64  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,3)[4], 52  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,3)[5], 40  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,3)[6], 68  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,3)[7], 56  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,MATCHED,3)[8], 44  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,0)[0], 39  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,0)[1], 51  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,0)[2], 63  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,0)[3], 43  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,0)[4], 55  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,0)[5], 67  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,0)[6], 47  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,0)[7], 59  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,0)[8], 71  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,1)[0], 63  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,1)[1], 67  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,1)[2], 71  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,1)[3], 51  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,1)[4], 55  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,1)[5], 59  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,1)[6], 39  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,1)[7], 43  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,1)[8], 47  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,2)[0], 71  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,2)[1], 59  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,2)[2], 47  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,2)[3], 67  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,2)[4], 55  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,2)[5], 43  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,2)[6], 63  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,2)[7], 51  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,2)[8], 39  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,3)[0], 47  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,3)[1], 43  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,3)[2], 39  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,3)[3], 59  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,3)[4], 55  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,3)[5], 51  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,3)[6], 71  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,3)[7], 67  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,MATCHED,3)[8], 63  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,0)[0], 72  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,0)[1], 84  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,0)[2], 96  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,0)[3], 76  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,0)[4], 88  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,0)[5], 100 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,0)[6], 80  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,0)[7], 92  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,0)[8], 104 );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,1)[0], 96  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,1)[1], 100 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,1)[2], 104 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,1)[3], 84  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,1)[4], 88  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,1)[5], 92  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,1)[6], 72  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,1)[7], 76  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,1)[8], 80  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,2)[0], 104 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,2)[1], 92  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,2)[2], 80  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,2)[3], 100 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,2)[4], 88  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,2)[5], 76  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,2)[6], 96  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,2)[7], 84  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,2)[8], 72  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,3)[0], 80  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,3)[1], 76  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,3)[2], 72  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,3)[3], 92  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,3)[4], 88  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,3)[5], 84  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,3)[6], 104 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,3)[7], 100 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,MATCHED,3)[8], 96  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,0)[0], 75  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,0)[1], 79  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,0)[2], 83  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,0)[3], 87  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,0)[4], 91  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,0)[5], 95  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,0)[6], 99  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,0)[7], 103 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,0)[8], 107 );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,1)[0], 83  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,1)[1], 95  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,1)[2], 107 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,1)[3], 79  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,1)[4], 91  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,1)[5], 103 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,1)[6], 75  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,1)[7], 87  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,1)[8], 99  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,2)[0], 107 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,2)[1], 103 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,2)[2], 99  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,2)[3], 95  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,2)[4], 91  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,2)[5], 87  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,2)[6], 83  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,2)[7], 79  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,2)[8], 75  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,3)[0], 99  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,3)[1], 87  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,3)[2], 75  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,3)[3], 103 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,3)[4], 91  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,3)[5], 79  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,3)[6], 107 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,3)[7], 95  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,MATCHED,3)[8], 83  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,0)[0], 0   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,0)[1], 4   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,0)[2], 8   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,0)[3], 12  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,0)[4], 16  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,0)[5], 20  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,0)[6], 24  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,0)[7], 28  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,0)[8], 32  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,1)[0], 24  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,1)[1], 12  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,1)[2], 0   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,1)[3], 28  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,1)[4], 16  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,1)[5], 4   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,1)[6], 32  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,1)[7], 20  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,1)[8], 8   );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,2)[0], 32  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,2)[1], 28  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,2)[2], 24  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,2)[3], 20  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,2)[4], 16  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,2)[5], 12  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,2)[6], 8   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,2)[7], 4   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,2)[8], 0   );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,3)[0], 8   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,3)[1], 20  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,3)[2], 32  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,3)[3], 4   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,3)[4], 16  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,3)[5], 28  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,3)[6], 0   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,3)[7], 12  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_NEG,INVERTED,3)[8], 24  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,0)[0], 3   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,0)[1], 15  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,0)[2], 27  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,0)[3], 7   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,0)[4], 19  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,0)[5], 31  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,0)[6], 11  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,0)[7], 23  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,0)[8], 35  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,1)[0], 11  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,1)[1], 7   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,1)[2], 3   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,1)[3], 23  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,1)[4], 19  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,1)[5], 15  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,1)[6], 35  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,1)[7], 31  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,1)[8], 27  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,2)[0], 35  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,2)[1], 23  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,2)[2], 11  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,2)[3], 31  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,2)[4], 19  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,2)[5], 7   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,2)[6], 27  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,2)[7], 15  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,2)[8], 3   );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,3)[0], 27  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,3)[1], 31  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,3)[2], 35  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,3)[3], 15  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,3)[4], 19  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,3)[5], 23  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,3)[6], 3   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,3)[7], 7   );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(KSI_POS,INVERTED,3)[8], 11  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[0], 36  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[1], 48  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[2], 60  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[3], 40  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[4], 52  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[5], 64  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[6], 44  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[7], 56  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[8], 68  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,1)[0], 44  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,1)[1], 40  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,1)[2], 36  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,1)[3], 56  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,1)[4], 52  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,1)[5], 48  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,1)[6], 68  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,1)[7], 64  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,1)[8], 60  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,2)[0], 68  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,2)[1], 56  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,2)[2], 44  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,2)[3], 64  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,2)[4], 52  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,2)[5], 40  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,2)[6], 60  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,2)[7], 48  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,2)[8], 36  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,3)[0], 60  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,3)[1], 64  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,3)[2], 68  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,3)[3], 48  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,3)[4], 52  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,3)[5], 56  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,3)[6], 36  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,3)[7], 40  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,3)[8], 44  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[0], 36  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[1], 48  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[2], 60  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[3], 40  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[4], 52  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[5], 64  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[6], 44  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[7], 56  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_NEG,INVERTED,0)[8], 68  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,0)[0], 39  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,0)[1], 43  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,0)[2], 47  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,0)[3], 51  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,0)[4], 55  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,0)[5], 59  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,0)[6], 63  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,0)[7], 67  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,0)[8], 71  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,1)[0], 63  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,1)[1], 51  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,1)[2], 39  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,1)[3], 67  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,1)[4], 55  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,1)[5], 43  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,1)[6], 71  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,1)[7], 59  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,1)[8], 47  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,2)[0], 71  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,2)[1], 67  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,2)[2], 63  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,2)[3], 59  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,2)[4], 55  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,2)[5], 51  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,2)[6], 47  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,2)[7], 43  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,2)[8], 39  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,3)[0], 47  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,3)[1], 59  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,3)[2], 71  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,3)[3], 43  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,3)[4], 55  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,3)[5], 67  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,3)[6], 39  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,3)[7], 51  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ETA_POS,INVERTED,3)[8], 63  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,0)[0], 72  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,0)[1], 76  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,0)[2], 80  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,0)[3], 84  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,0)[4], 88  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,0)[5], 92  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,0)[6], 96  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,0)[7], 100 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,0)[8], 104 );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,1)[0], 96  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,1)[1], 84  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,1)[2], 72  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,1)[3], 100 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,1)[4], 88  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,1)[5], 76  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,1)[6], 104 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,1)[7], 92  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,1)[8], 80  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,2)[0], 104 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,2)[1], 100 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,2)[2], 96  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,2)[3], 92  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,2)[4], 88  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,2)[5], 84  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,2)[6], 80  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,2)[7], 76  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,2)[8], 72  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,3)[0], 80  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,3)[1], 92  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,3)[2], 104 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,3)[3], 76  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,3)[4], 88  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,3)[5], 100 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,3)[6], 72  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,3)[7], 84  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_NEG,INVERTED,3)[8], 96  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,0)[0], 75  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,0)[1], 87  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,0)[2], 99  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,0)[3], 79  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,0)[4], 91  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,0)[5], 103 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,0)[6], 83  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,0)[7], 95  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,0)[8], 107 );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,1)[0], 83  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,1)[1], 79  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,1)[2], 75  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,1)[3], 95  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,1)[4], 91  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,1)[5], 87  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,1)[6], 107 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,1)[7], 103 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,1)[8], 99  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,2)[0], 107 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,2)[1], 95  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,2)[2], 83  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,2)[3], 103 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,2)[4], 91  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,2)[5], 79  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,2)[6], 99  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,2)[7], 87  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,2)[8], 75  );

  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,3)[0], 99  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,3)[1], 103 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,3)[2], 107 );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,3)[3], 87  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,3)[4], 91  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,3)[5], 95  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,3)[6], 75  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,3)[7], 79  );
  BOOST_CHECK_EQUAL( sf->face_flx_pts(ZTA_POS,INVERTED,3)[8], 83  );
}

////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE( finalize_mpi )
{
#ifdef test_is_mpi
  PE::Comm::instance().finalize();
#endif
}

///////////////////////////////////////////////////////////////////////////////



BOOST_AUTO_TEST_SUITE_END()

////////////////////////////////////////////////////////////////////////////////
