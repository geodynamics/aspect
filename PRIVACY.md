Privacy Policy
==============

We as the developers of the ASPECT project collect basic information about the
use of ASPECT using Google Analytics. This information is of great help for us
to further develop and improve ASPECT, e.g. to make decisions about versions in
use and backwards compatibility. It also allows us to justify ASPECT's impact
to our funding agencies, e.g. by estimating the size of our user community,
which is notoriously difficult for open-source community projects like ASPECT.
We are aware of the justified concerns about user privacy and therefore adopt
the following policies: 

We collect **only** the following information, which is also displayed on
screen when ASPECT starts:

* ASPECT version number (excluding git branch or commit hash),
* Version number of the main dependencies (deal.II, trilinos, p4est),
* Limited compile-time configurations (build type and vectorization level),
* Number of MPI processes used.

All source code collecting this information can be found in 'source/main.cc' in
the function 'track_usage()'.

Google Analytics also automatically collects an anonymized version of your IP
for geolocational purposes. This IP can not be traced back to you, but does
allow the identification of your country and potentially your city. While this
is not our intent, we cannot prevent the collection of this data.


Participation
-------------

You may opt out of this tracking by setting the cmake variable
'ASPECT_TRACK_USAGE_DATA=OFF' during configuration, which will disable all
tracking code. We ask that you only do so if you have serious concerns about
the collected data (in which case please also send us an email to discuss your
concerns, we are open to improvements). Providing this data is a valuable
contribution that supports ASPECT's development.


Policy
------

1. We respect the privacy of our users and will never collect information
   that can be used to identify individual users or that can be used to
   reconstruct the intent of your use of ASPECT.

2. We collect information only where this is necessary for the above specified
   purpose.

3. We will only use this information for the purpose for which it was 
   collected.
