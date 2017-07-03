#include "math/Legendre.h"
#include "BoostUnitTest.h"

BOOST_AUTO_TEST_CASE(LegendreRoots) {
    std::vector< std::vector <double> > legendreRootsKnown{
        //order 1
        {
            0.0
        },
        //order 2
        {
            -0.57735026918962584,
            0.57735026918962562
        },
        //order 3
        {
            -0.7745966692414834,
            0.0,
            0.77459666924148352
        },
        //order 4
        {
            -0.86113631159405335,
            -0.3399810435848562,
            0.33998104358485665,
            0.86113631159405257
        },
        //order 5
        {
            -0.9061798459386643,
            -0.53846931010568333,
            0.0,
            0.53846931010568311,
            0.90617984593866352
        },
        //order 6
        {
            -0.93246951420315138,
            -0.66120938646626437,
            -0.23861918608319704,
            0.23861918608319688,
            0.66120938646626337,
            0.93246951420315194
        },
        //order 7
        {
            -0.94910791234275782,
            -0.74153118559939579,
            -0.40584515137739685,
            0.0,
            0.40584515137739685,
            0.74153118559939546,
            0.94910791234275971
        },
        //order 8
        {
            -0.9602898564975284,
            -0.79666647741363394,
            -0.52553240991632688,
            -0.1834346424956497,
            0.18343464249564989,
            0.52553240991632955,
            0.79666647741362762,
            0.96028985649753373
        },
        //order 9
        {
            -0.9681602395076313,
            -0.83603110732662467,
            -0.6133714327005938,
            -0.32425342340380914,
            0.0,
            0.32425342340380908,
            0.61337143270059147,
            0.83603110732663521,
            0.96816023950762675
        }
    };

  int maxOrder = 9;

  for (int order = 1; order <= maxOrder; order++) {
    std::vector<double> roots = NuTo::Math::Polynomial::LegendreRoots(order);
    for (int i=0; i<order; i++) {
        BOOST_CHECK_CLOSE(roots[i],legendreRootsKnown[order-1][i],1e-10);
    }
  }
}

BOOST_AUTO_TEST_CASE(LegendreDerivRoots) {
    std::vector< std::vector <double> > legendreDerivRootsKnown{
        //order 1
            {},
        //order 2
            {0.0},
        //order 3
            {
                -0.44721359549995798,
                0.44721359549995798
            },
        //order 4
            {
                -0.65465367070797731,
                3.2381504884900401e-17,
                0.65465367070797731
            },
        //order 5
            {
                -0.76505532392946518,
                -0.28523151648064488,
                0.28523151648064504,
                0.76505532392946551
            },
        //order 6
            {
                -0.83022389627856719,
                -0.46884879347071395,
                0.0,
                0.46884879347071434,
                0.83022389627856663
            },
        //order 7
            {
                -0.87174014850960502,
                -0.59170018143314251,
                -0.2092992179024788,
                0.20929921790247896,
                0.5917001814331424,
                0.87174014850960502
            },
        //order 8
            {
                -0.89975799541145662,
                -0.67718627951074062,
                -0.36311746382617754,
                -2.214870152363245e-18,
                0.36311746382617821,
                0.67718627951073862,
                0.89975799541146051
            },
        //order 9
            {
                -0.91953390816645975,
                -0.7387738651055038,
                -0.47792494981044498,
                -0.16527895766638714,
                0.16527895766638703,
                0.47792494981044326,
                0.73877386510550713,
                0.91953390816645708
            }
    };

  int maxOrder = 9;

  for (int order = 2; order <= maxOrder; order++) {
    std::vector<double> roots = NuTo::Math::Polynomial::LegendreDerivRoots(order);
    for (int i=0; i<(order-1); i++) {
        if ( (roots[i] > 1e-10 ) && (legendreDerivRootsKnown[order-1][i] > 1e-10 ))
            BOOST_CHECK_CLOSE(roots[i],legendreDerivRootsKnown[order-1][i],1e-10);
    }
  }
}
