print "** Lua ephemeris validation test **"

-- Convert degrees, minutes, seconds to decimal
function dms(deg, min, sec)
	return deg + min/60 + sec/3600
	end


planet = Sun
Tsnapshot = {d=1, m=8, y=2018}

print ("Geo longitude = ", dms(128, 27, 56.86))
print ("Geo latitude  = ", dms(  0,  0, -6.87))
print  "Rayon vector    = 1.01500309 ua"

