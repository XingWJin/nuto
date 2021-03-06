#pragma once

#include <vector>
#include <eigen3/Eigen/Core>

namespace NuTo
{
namespace Visualize
{
//! @brief ... point for visualization
class Point
{
public:
    //! @brief ... constructor
    //! @param coordinates ... point coordinates
    //! @param numData ... number of different point data fields
    Point(Eigen::Vector3d coordinates, int numData);

    //! @brief ... set data
    //! @param dataIndex ... data index
    //! @param data ... data
    void SetData(int dataIndex, Eigen::VectorXd data);

    //! @brief ... get point coordinates
    //! @return ... reference to point coordinates
    const Eigen::Vector3d& GetCoordinates() const;

    //! @brief ... get data
    //! @param dataIndex ... data index
    //! @return ... reference to point data
    const Eigen::VectorXd& GetData(int dataIndex) const;

private:
    //! @brief ... point coordinates
    Eigen::Vector3d mCoordinates;

    //! @brief ... vector of point data
    std::vector<Eigen::VectorXd> mData;
};
} // namespace Visualize
} // namespace NuTo
