# The pass an optical scene was imaged on, from its satellite's trajectory.

using Dates: DateTime
using AutoriftItsLiveProduction: landsat_flight_direction, s2_flight_direction

# An `_ANG.txt` ephemeris block of four points whose ECEF z changes by `dz` per second.
ang(dz) = """
GROUP = EPHEMERIS
  NUMBER_OF_POINTS = 4
  EPHEMERIS_TIME = (0.0, 1.0, 2.0, 3.0)
  EPHEMERIS_ECEF_X = (1.0, 2.0, 3.0, 4.0)
  EPHEMERIS_ECEF_Y = (1.0, 2.0, 3.0, 4.0)
  EPHEMERIS_ECEF_Z = ($(6e6), $(6e6 + dz), $(6e6 + 2dz), $(6e6 + 3dz))
END_GROUP = EPHEMERIS
"""

gps(t, vz) = "<GPS_Point><POSITION_VALUES unit=\"mm\">1 2 3</POSITION_VALUES>" *
             "<VELOCITY_VALUES unit=\"mm/s\">10 20 $vz</VELOCITY_VALUES><GPS_TIME>$t</GPS_TIME></GPS_Point>"

@testset "a Landsat scene's pass is the sign of its satellite's northward motion" begin
    @test landsat_flight_direction(ang(-3600.0)) == "descending"
    @test landsat_flight_direction(ang(3600.0)) == "ascending"
end

@testset "a Sentinel-2 tile's pass is its own, not its datastrip's start" begin
    # A datastrip that starts northbound before the pole and is southbound by the tile's time.
    head = "<Datastrip><GPS_Points_List>" * gps("2020-07-01T20:40:06.0", 256608) *
           gps("2020-07-01T20:41:06.0", -203656) * gps("2020-07-01T20:48:30.0", -3468420) *
           "</GPS_Points_List>"
    @test s2_flight_direction(head, DateTime(2020, 7, 1, 20, 48, 30)) == "descending"
    @test s2_flight_direction(head, DateTime(2020, 7, 1, 20, 40, 10)) == "ascending"
    @test_throws "no GPS fix within a minute" s2_flight_direction(head, DateTime(2020, 7, 1, 21, 0, 0))
end
