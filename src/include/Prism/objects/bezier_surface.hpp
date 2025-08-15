#ifndef PRISM_BEZIER_SURFACE_HPP_
#define PRISM_BEZIER_SURFACE_HPP_

#include "prism_export.h"

#include "Prism/core/material.hpp"
#include "Prism/core/point.hpp"
#include "Prism/core/ray.hpp"
#include "Prism/core/vector.hpp"
#include "Prism/objects/objects.hpp"
#include "Prism/objects/triangle.hpp"

#include <memory>
#include <vector>

namespace Prism {

/**
 * @class BezierSurface
 * @brief Represents a Bezier surface in 3D space defined by control points arranged in a grid.
 * This class inherits from Object and implements the hit method to check for ray-surface
 * intersections using adaptive subdivision and numerical methods.
 */
class BezierSurface : public Object {
  public:
    /**
     * @brief Constructs a BezierSurface given a 2D grid of control points and a material.
     * @param control_points A 2D vector of Point3 representing the control points grid.
     *                      control_points[i][j] represents the control point at row i, column j.
     * @param material A shared pointer to a Material object associated with the surface.
     * @param u_segments Number of segments in the u direction for rendering.
     * @param v_segments Number of segments in the v direction for rendering.
     */
    BezierSurface(const std::vector<std::vector<Point3>>& control_points, 
                  std::shared_ptr<Material> material, int u_segments, int v_segments);


    /**
     * @brief Checks if a ray intersects with the Bezier surface.
     * @param ray The Ray to test for intersection.
     * @param t_min The minimum distance for a valid hit.
     * @param t_max The maximum distance for a valid hit.
     * @param rec The HitRecord to be filled upon a collision.
     * @return True if the ray intersects the surface within the specified distance range, false
     * otherwise.
     */
    virtual bool hit(const Ray& ray, double t_min, double t_max, HitRecord& rec) const override;

    /**
     * @brief Gets the number of control curves in the u direction.
     * @return The number of control points in the u direction.
     */
    int getUCount() const { return u_count; }

    /**
     * @brief Gets the number of control curves in the v direction.
     * @return The number of control points in the v direction.
     */
    int getVCount() const { return v_count; }

    /**
     * @brief Evaluates the Bezier surface at parameter coordinates (u, v).
     * @param u Parameter in the u direction (0.0 to 1.0).
     * @param v Parameter in the v direction (0.0 to 1.0).
     * @return The point on the surface at the given parameters.
     */
    Point3 evaluateAt(double u, double v) const;

    /**
     * @brief Computes the normal vector at parameter coordinates (u, v).
     * @param u Parameter in the u direction (0.0 to 1.0).
     * @param v Parameter in the v direction (0.0 to 1.0).
     * @return The normalized normal vector at the given parameters.
     */
    Vector3 normalAt(double u, double v) const;

  private:
    std::vector<std::vector<Point3>> control_points; ///< 2D grid of control points
    int u_count; ///< Number of control points in u direction
    int v_count; ///< Number of control points in v direction
    std::shared_ptr<Material> material; ///< Material properties of the surface
    
    // Mesh representation for fast ray intersection
    std::vector<MeshTriangle> triangulated_mesh; ///< Tessellated triangles representing the surface
    
    // Bounding box for early ray rejection
    Point3 bbox_min; ///< Minimum corner of bounding box
    Point3 bbox_max; ///< Maximum corner of bounding box
    
    // Tessellation parameters
    int u_segments; ///< Number of segments in u direction for tessellation
    int v_segments; ///< Number of segments in v direction for tessellation

    /**
     * @brief Computes Bernstein polynomial basis function.
     * @param i Index of the basis function.
     * @param n Degree of the polynomial.
     * @param t Parameter value.
     * @return The value of the Bernstein basis function.
     */
    double bernstein(int i, int n, double t) const;

    /**
     * @brief Computes binomial coefficient C(n, k).
     * @param n Upper value.
     * @param k Lower value.
     * @return The binomial coefficient.
     */
    int binomial(int n, int k) const;

    /**
     * @brief Tessellates the Bezier surface into triangles for fast intersection.
     * This method evaluates the surface at regular intervals and creates triangles.
     */
    void tessellateSurface();

    /**
     * @brief Computes the bounding box of the surface.
     * Updates bbox_min and bbox_max based on the tessellated mesh.
     */
    void computeBoundingBox();

    /**
     * @brief Tests if a ray intersects the bounding box.
     * @param ray The ray to test.
     * @param t_min Minimum distance.
     * @param t_max Maximum distance.
     * @return True if ray intersects bounding box.
     */
    bool rayIntersectsBoundingBox(const Ray& ray, double t_min, double t_max) const;
};

} // namespace Prism

#endif // PRISM_BEZIER_SURFACE_HPP_
