
#include <gtest/gtest.h>

#include <motor_mount.h>

class MotorMountTest : public testing::Test {
protected:
    Motor motor; // member variable visible to all TEST_F cases
    MotorMount motor_mount;

    MotorMountTest()
        : motor(12.0, 2.4, 0.0), motor_mount(0.0, 0.0, 10.0, 0.0, Eigen::Vector3d(-0.5, 0.0, 0.0), motor)
    {}   // initialize here

    ~MotorMountTest() override = default;

    void SetUp() override {
        double t = 0.0;
        motor_mount.Ignite(t);
    }

    void TearDown() override {
        //
    }
};

TEST_F(MotorMountTest, TestPositionVector) {
    double t = 0.0;
    motor_mount.Command(0.0, 0.0, t);
    Eigen::Vector3d position = motor_mount.GetPositionVector();

    EXPECT_TRUE(std::abs(position.x() + 0.5) <= 1e-3);
    EXPECT_TRUE(std::abs(position.y()) <= 1e-3);
    EXPECT_TRUE(std::abs(position.z()) <= 1e-3);
}

TEST_F(MotorMountTest, TestThrustMoment_Zero) {
    double t = 0.0;
    motor_mount.Command(0.0, 0.0, t);
    Eigen::Vector3d thrust = motor_mount.GetThrustVector(t);
    Eigen::Vector3d position = motor_mount.GetPositionVector();
    Eigen::Vector3d moment = position.cross(thrust);

    EXPECT_TRUE(std::abs(moment.x()) <= 1e-3);
    EXPECT_TRUE(std::abs(moment.y()) <= 1e-3);
    EXPECT_TRUE(std::abs(moment.z()) <= 1e-3);
}

TEST_F(MotorMountTest, TestThrustMoment_X) {
    double t = 0.0;
    motor_mount.Command(90.0, 0.0, t);
    Eigen::Vector3d thrust = motor_mount.GetThrustVector(t);
    Eigen::Vector3d position = motor_mount.GetPositionVector();
    Eigen::Vector3d moment = position.cross(thrust);

    EXPECT_TRUE(std::abs(moment.x()) <= 1e-3);
    EXPECT_TRUE(std::abs(moment.y()) <= 1e-3);
    EXPECT_TRUE(std::abs(moment.z() + 6.0) <= 1e-3);
}

TEST_F(MotorMountTest, TestThrustMoment_MinusX) {
    double t = 0.0;
    motor_mount.Command(-90.0, 0.0, t);
    Eigen::Vector3d thrust = motor_mount.GetThrustVector(t);
    Eigen::Vector3d position = motor_mount.GetPositionVector();
    Eigen::Vector3d moment = position.cross(thrust);

    EXPECT_TRUE(std::abs(moment.x()) <= 1e-3);
    EXPECT_TRUE(std::abs(moment.y()) <= 1e-3);
    EXPECT_TRUE(std::abs(moment.z() - 6.0) <= 1e-3);
}

TEST_F(MotorMountTest, TestThrustMoment_Y) {
    double t = 0.0;
    motor_mount.Command(0.0, 90.0, t);
    Eigen::Vector3d thrust = motor_mount.GetThrustVector(t);
    Eigen::Vector3d position = motor_mount.GetPositionVector();
    Eigen::Vector3d moment = position.cross(thrust);

    EXPECT_TRUE(std::abs(moment.x()) <= 1e-3);
    EXPECT_TRUE(std::abs(moment.y() - 6.0) <= 1e-3);
    EXPECT_TRUE(std::abs(moment.z()) <= 1e-3);
}

TEST_F(MotorMountTest, TestThrustMoment_MinusY) {
    double t = 0.0;
    motor_mount.Command(0.0, -90.0, t);
    Eigen::Vector3d thrust = motor_mount.GetThrustVector(t);
    Eigen::Vector3d position = motor_mount.GetPositionVector();
    Eigen::Vector3d moment = position.cross(thrust);

    EXPECT_TRUE(std::abs(moment.x()) <= 1e-3);
    EXPECT_TRUE(std::abs(moment.y() + 6.0) <= 1e-3);
    EXPECT_TRUE(std::abs(moment.z()) <= 1e-3);
}

TEST_F(MotorMountTest, TestThrustVector_Straight) {
    double t = 0.0;
    motor_mount.Command(0.0, 0.0, t);
    Eigen::Vector3d thrust = motor_mount.GetThrustVector(t);

    EXPECT_TRUE(std::abs(thrust.x() - 12.0) <= 1e-3);
    EXPECT_TRUE(std::abs(thrust.y()) <= 1e-3);
    EXPECT_TRUE(std::abs(thrust.z()) <= 1e-3);
}

TEST_F(MotorMountTest, TestThrustVector_FullX) {
    double t = 0.0;
    motor_mount.Command(90.0, 0.0, t);
    Eigen::Vector3d thrust = motor_mount.GetThrustVector(t);

    EXPECT_TRUE(std::abs(thrust.x()) <= 1e-3);
    EXPECT_TRUE(std::abs(thrust.y() - 12.0) <= 1e-3);
    EXPECT_TRUE(std::abs(thrust.z()) <= 1e-3);
}

TEST_F(MotorMountTest, TestThrustVector_FullY) {
    double t = 0.0;
    motor_mount.Command(0.0, 90.0, t);
    Eigen::Vector3d thrust = motor_mount.GetThrustVector(t);

    EXPECT_TRUE(std::abs(thrust.x()) <= 1e-3);
    EXPECT_TRUE(std::abs(thrust.y()) <= 1e-3);
    EXPECT_TRUE(std::abs(thrust.z() - 12.0) <= 1e-3);
}

TEST_F(MotorMountTest, TestBurnout) {
    double t0 = 0.0, tf = 2.4, tl = 2.41;

    Eigen::Vector3d thrust;

    thrust = motor_mount.GetThrustVector(t0);
    EXPECT_DOUBLE_EQ(thrust.x(), 12.0);
    EXPECT_DOUBLE_EQ(thrust.y(), 0.0);
    EXPECT_DOUBLE_EQ(thrust.z(), 0.0);

    thrust = motor_mount.GetThrustVector(tf);
    EXPECT_DOUBLE_EQ(thrust.x(), 12.0);
    EXPECT_DOUBLE_EQ(thrust.y(), 0.0);
    EXPECT_DOUBLE_EQ(thrust.z(), 0.0);

    thrust = motor_mount.GetThrustVector(tl);
    EXPECT_DOUBLE_EQ(thrust.x(), 0.0);
    EXPECT_DOUBLE_EQ(thrust.y(), 0.0);
    EXPECT_DOUBLE_EQ(thrust.z(), 0.0);
}
