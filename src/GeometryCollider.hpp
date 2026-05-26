#pragma once
#include <psDomain.hpp>

template <int DIM>
class GeometryCollider {
    private:
        RTCDevice device;
        RTCScene scene;

        struct point_4f_t { float xx, yy, zz, radius; };
        struct normal_vec_3f_t { float xx, yy, zz; };

    public:
        // Constructor
        GeometryCollider() {
            device = rtcNewDevice("hugepages=1");
            scene = rtcNewScene(device);
            rtcSetSceneFlags(scene, RTC_SCENE_FLAG_NONE);
            rtcSetSceneBuildQuality(scene, RTC_BUILD_QUALITY_HIGH);
        }

        // Destructor
        ~GeometryCollider() {
            rtcReleaseScene(scene);
            rtcReleaseDevice(device);
        }

        // Loade ViennaPS DiskMesh into Embree Scene
        void buildSurface(viennacore::SmartPointer<viennals::Mesh<double>> mesh, double gridDelta) {
            auto &points = mesh->nodes;
            auto &normals = *mesh->getCellData().getVectorData("Normals");
            auto numPoints = points.size();

            double diskRadii = gridDelta * rayInternal::DiskFactor<DIM>;

            RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_ORIENTED_DISC_POINT);

            // Initialize Embree geometry buffers for vertices and normals
            point_4f_t *pointBuf = (point_4f_t *)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT4, sizeof(point_4f_t), numPoints);
            normal_vec_3f_t *normalBuf = (normal_vec_3f_t *)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_NORMAL, 0, RTC_FORMAT_FLOAT3, sizeof(normal_vec_3f_t), numPoints);
            
            // Fill the geometry buffers with mesh data
            for (size_t i = 0; i < numPoints; ++i) {
                pointBuf[i] = {(float)points[i][0], (float)points[i][1], (float)points[i][2], (float)diskRadii};
                normalBuf[i] = {(float)normals[i][0], (float)normals[i][1], 0.f};
            }

            rtcCommitGeometry(geom);
            rtcAttachGeometry(scene, geom);
            rtcReleaseGeometry(geom);
            rtcCommitScene(scene);
        }

        // Shoot a ray; return NodeID if hit, -1 if miss
        int checkCollision(const std::array<double, 3>& origin, const std::array<double, 3>& direction, double max_distance) {
            RTCRayHit rayHit;
            rayHit.ray.org_x = (float)origin[0]; 
            rayHit.ray.org_y = (float)origin[1]; 
            rayHit.ray.org_z = 0.f;
            
            rayHit.ray.dir_x = (float)direction[0]; 
            rayHit.ray.dir_y = (float)direction[1]; 
            rayHit.ray.dir_z = 0.f;
            
            rayHit.ray.tnear = 0.0f;
            rayHit.ray.tfar = (float)max_distance;
            rayHit.ray.mask = -1;
            rayHit.ray.flags = 0;
            rayHit.hit.geomID = RTC_INVALID_GEOMETRY_ID;

            rtcIntersect1(scene, &rayHit);

            if (rayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
                return rayHit.hit.primID;
            }
            return -1;
        }
};