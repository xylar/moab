#include <math.h>
#include <vector>
#include <stdexcept>
#include <assert.h>

class geomObject
{
    public:
        virtual ~geomObject() {}
        virtual void project_points2geom ( int dim, double* oldcoords, double* newcoords, double* derivs ) const = 0;
        virtual double compute_projecterror ( int dim, double* oldcoords ) const
        {
            double* newcoords = new double[dim]();
            project_points2geom ( dim, oldcoords, newcoords, NULL );
            double ans = 0;

            for ( int i = 0; i < dim; ++i ) { ans += ( newcoords[i] - oldcoords[i] ) * ( newcoords[i] - oldcoords[i] ); }

            delete[] newcoords;
            return sqrt ( ans );
        }
        virtual void compute_projecterror ( int dim, int nverts, double* oldcoords, double& l1err, double& l2err, double& linferr ) const
        {
            l1err = l2err = linferr = 0;

            for ( int i = 0; i < nverts; ++i )
            {
                double err = compute_projecterror ( dim, oldcoords + i * dim );
                l1err += err; l2err += err * err; linferr = std::max ( linferr, err );
            }

            l1err /= nverts; l2err = sqrt ( l2err / nverts );
        }
        inline double Twonorm ( int dim, double* vec ) const
        {
            double ans = 0;

            for ( int i = 0; i < dim; ++i ) { ans += vec[i] * vec[i]; }

            return sqrt ( ans );
        }
};

class sphere : public geomObject
{
    public:
        sphere ( double x = 0, double y = 0, double z = 0, double r = 1 ) : centerx ( x ), centery ( y ), centerz ( z ), radius ( r ) {assert ( r > 0 );}
        virtual ~sphere() {}
        void project_points2geom ( int dim, double* oldcoords, double* newcoords, double* derivs ) const
        {
            if ( oldcoords == NULL || newcoords == NULL ) { throw std::invalid_argument ( "NULL pointer" ); }

            std::vector<double> direction;
            direction.push_back ( oldcoords[0] - centerx ); direction.push_back ( oldcoords[1] - centery );

            if ( dim == 3 ) { direction.push_back ( oldcoords[2] - centerz ); }

            double len = geomObject::Twonorm ( dim, & ( direction[0] ) ); assert ( len > 0 );

            for ( int i = 0; i < dim; ++i ) { direction[i] /= len; }

            newcoords[0] = centerx + direction[0] * radius; newcoords[1] = centery + direction[1] * radius;

            if ( dim == 2 )
            {
                if ( derivs )
                {
                    derivs[0] = -direction[1]; derivs[1] = direction[0];
                }

                return;
            }
            else if ( dim == 3 )
            {
                newcoords[2] = centerz + direction[2] * radius;

                if ( derivs )
                {
                    derivs[0] = direction[0]; derivs[1] = direction[1]; derivs[2] = direction[2];
                }
            }
            else
            {
                throw std::invalid_argument ( "dim must be 2 or 3" );
            }
        }
    private:
        double centerx, centery, centerz, radius;
};

class torus : public geomObject
{
    public:
        torus ( double x = 0, double y = 0, double z = 0, double aa = 0.3, double cc = 1.0 ) : centerx ( x ), centery ( y ), centerz ( z ), a ( aa ), c ( cc ) {assert ( aa > 0 && cc > 0 );}
        virtual ~torus() {}
        void project_points2geom ( int dim, double* oldcoords, double* newcoords, double* derivs ) const
        {
            assert ( dim == 3 && oldcoords && newcoords );
            double transfer[3] = {oldcoords[0] - centerx, oldcoords[1] - centery, oldcoords[2] - centerz};
            double twodnrm = geomObject::Twonorm ( 2, transfer );
            double tubecenter[3] = {c* transfer[0] / twodnrm + centerx, c* transfer[1] / twodnrm + centery, centerz};
            double direction[3] = {oldcoords[0] - tubecenter[0], oldcoords[1] - tubecenter[1], oldcoords[2] - tubecenter[2]};
            double len = geomObject::Twonorm ( 3, direction ); assert ( len > 0 );
            direction[0] /= len; direction[1] /= len; direction[2] /= len;

            if ( derivs )
            {
                derivs[0] = direction[0]; derivs[1] = direction[1]; derivs[2] = direction[2];
            }

            newcoords[0] = tubecenter[0] + a * direction[0]; newcoords[1] = tubecenter[1] + a * direction[1]; newcoords[2] = tubecenter[2] + a * direction[2];
        }
    private:
        double centerx, centery, centerz, a, c;
};

