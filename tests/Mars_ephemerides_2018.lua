print "** Lua ephemeris validation test **"

-- Convert degrees, minutes, seconds to decimal
function dms(deg, min, sec)
	return deg + min/60 + sec/3600
	end

planet = Mars
Tsnapshot = {d=28, m=1, y=2018}

print ("Helio longitude = ", dms(206, 16, 18.51))
print ("Helio latitude  = ", dms(  0, 43, 45.92))
print  "Rayon vector    = 1.60662322 ua"

