# Ball Flight

Ballflight enables the modeling of a pitched or hit ball. It is useful in determining the results of the pitch or the
hit, and also includes plotting functions that are helpful in visualizing the metric outputs of the ballflight model.


## Install
```
$ pip install ballflight
```

## Command Line Instructions
For a full list of commands type
```
$ python ballflight --help
```
To get the results of a particular pitch:
```
$ python ballflight --releaseSpeed 88.4 --spinRate 2589 --breakSpin 2589 --spinDirection 12:33
```
To plot the results of a particular pitch, add a `-p` option to the end of your command line:
```
$ python ballflight --releaseSpeed 88.4 --spinRate 2589 --breakSpin 2589 --spinDirection 12:33 -p
```

## Usage Instructions

In order to use Ballflight, first create a `PitchOptions` or `SwingOptions` object with your input information.

For example:
```
pitch_options = PitchOptions(play_type='FAST_PITCH_SOFTBALL', release_speed=60.0, release_spin_rate=900,
                             release_spin_axis=[1, 0, 0], release_lauch_angle=5.6, release_heading_angle=2.0, pitch=True)
```

Then you call a `PitchResults` object as follows:
```
pitch_result = PitchResult.throw_ball(pitch_options)
```





## Release History
- v0.1.1, October 19, 2018
  - __main__ method for ball flight with command line args.
- v0.1.0, August 9, 2018
  - Initial Release


## License

MIT (see `LICENSE` for more information).


## Thank You!

**Many** thanks to the Diamond Kinetics team, especially Minmin Zhang and Mike Ressler.


## Works Referenced

The following works were consulted in the creation of this model.
PRIMARY SOURCES CONSULTED:
1. Baseball-bat collisions and the resulting trajectories of spinning balls -- Robert G. Watts and Steven Baroni
2. [How to hit home runs: Optimum baseball bat swing parameters for maximum range trajectories -- Gregory S. Sawicki, Mont Hubbard, and William J. Stronge](https://pdfs.semanticscholar.org/3311/8293e2a5f7fadde1cbdc6c5a9a67b10706fa.pdf)
3. Spin of a batted baseball -- Alan M. Nathan, Jonas Cantakos, Russ Kesman, Biju Mathew, Wes Lukash
    http://baseball.physics.illinois.edu/ProcediaEngineering34Spin.pdf
4. Aerodynamic drag crisis and its possible effect on the flight of baseballs -- Cliff Frolich
5. List of moments of inertia -- https://en.wikipedia.org/wiki/List_of_moments_of_inertia
6. Personal Correspondence with Mont Hubbard
7. Report of the Committee Studying Home Run Rates in Major League Baseball -- Jim Albert, Jay Bartroff,
    Roger Blandford , Dan Brooks, Josh Derenski, Larry Goldstein, Anette (Peko) Hosoi, Gary Lorden, Alan Nathan,
    and Lloyd Smith

OTHER SOURCES:
- Dynamics of the baseball-ball collision -- Alan M. Nathan
- Comment on Source 2 (Robert Adair) and Reply to the Comment (Sawicki, Hubbard, Stronge)
- The effect of spin on the flight of batted baseballs - A. F. Rex
- The effect of spin on the flight of a baseball - Alan M. Nathan
- The physics of baseball - Robert Adair
- The effects of surface roughness and tunnel blockage on the flow past spheres -- Elmar Achenbach
- Determining aerodynamic properties of sports balls in situ -- Jeffrey Kensrud
- Observed drag crisis on a sphere in flowing He I and He II -- Michael R. Smith, David K. Hilton, Steven W. Van Sciver
- Baseball (ball) - https://en.wikipedia.org/wiki/Baseball_(ball)
- [Coefficient of Restitution: A Comparison of Major League and Little League Baseballs -- Cameron Wallace](http://cssf.usc.edu/History/2008/Projects/J1935.pdf)
- Softball - https://en.wikipedia.org/wiki/Softball
- NCAA Softball Rules - http://fs.ncaa.org/Docs/stats/Stats_Manuals/Softball/Softball_Rules.pdf

