print "** Lua ephemeris validation test **"

-- Convert degrees, minutes, seconds to decimal
function dms(deg, min, sec)
	return deg + min/60 + sec/3600
	end


planet = Jupiter
Tsnapshot = {d=5, m=3, y=2018}

print ("Helio longitude = ", dms(223,  7, 18.55))
print ("Helio latitude  = ", dms(  1,  5, 53.81))
print  "Rayon vector    = 5.4218701 ua"

