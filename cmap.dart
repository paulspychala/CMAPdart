library cmap;

/*
 * Causal Mapping Library written in Dart for use in VisiNets
 * Copyright (c) 2015 VisiNets Inc
 * 
 * Literature:
 * Causal mapping as a tool to mechanistically interpret phenomena in cell motility: Application to cortical oscillations in spreading cells
 * - http://onlinelibrary.wiley.com/doi/10.1002/cm.20143/abstract
 * In Silico Generation of Alternative Hypotheses Using Causal Mapping (CMAP)
 * - http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0005378
 * 
 * Home Info:
 * http://www.visinets.com
 */

import "dart:math";

class Species {
  // name of species
  String name;
  
  // values
  double previous_value = 0.0;
  double future_value = 0.0;
  
  // alpha value
  double alpha;
  num cmax;
  
  // parameters
  Map parameters;
  
  // lock species value from updating
  bool lock = false;
  
  // list of reactions
  List<Influence> reactionsChild;
  List<Influence> reactionsParent;
  
  // constructor
  Species(){   
   reactionsChild = new List<Influence>();
   reactionsParent = new List<Influence>();
  
   parameters = new Map();
   parameters['inherentActivity'] = 0.0;
  }
}

class Influence {
  // temp value for storing value of influence
  double value;
  
  // influences have several parameters:
  // - weight
  // - activation
  // - inhibition
  Map parameters;
  
  List<Species> to;
  List<Species> from;
  
  // constructor
  Influence(){
    to = new List<Species>();
    from = new List<Species>();
    
    parameters = new Map();
    parameters['weight'] = .5;
      
  }
}

class Network {
  // list of species in the network
  List<Species> nodes;
  
  // how many runs in the simulation
  num runLength;
  
  // constructor
  Network(){
    nodes = new List<Species>();
  }
}

class CMapSimulation {
  int currentTime = 0;

  final int MAX_SIZE = 1;

  bool calculateSpecies(Species p)
  {
    // x = a Σ w C(t-1)
    double x = 0.0;
    Influence l;
    for(int i = 0; i < p.reactionsParent.length; i++)
    {
      l = p.reactionsParent[i];
      l.value = l.parameters['weight'];
      for(var sp in l.from){
        l.value = l.value * sp.previous_value;
      }
      if(l.parameters['activation'] == true)
      {
        x += l.value;
      } else if(l.parameters['inhibition'] == true) {
        x -= l.value;
      }
    }
    x = x * p.alpha;

    // f = 1 - e^-x / 1 + e^-x
    double f = (1 - exp(x * -1)) / (1 + exp(x * -1));

    // OLD -> A = 1-C(t-1) if f > 0
    // NEW -> A = Cmax - C(t-1) if f > 0
    // A = C(t-1) if f <= 0
    double A;
    if(f > 0) {
      A = p.cmax - p.previous_value;
    } else {
      A = p.previous_value;
    }

    // C = C(t-1) + A(c,f)f(x)
    double C = p.previous_value + A * f;

    if(!p.lock) {
      p.future_value = C;
    }

    return true;
  }
  

  bool runAll(Network network){
    for(var s in network.nodes){
      s.future_value = s.parameters['inherentActivity'];
    }

    for(currentTime = 0; currentTime < network.runLength; currentTime++){
      // prepare run
      for(var s in network.nodes) {
        s.previous_value = s.future_value;
      }

      // calc newrun
      for(var s in network.nodes){
        calculateSpecies(s);
      }
    }
    
    return true;
  }

}
