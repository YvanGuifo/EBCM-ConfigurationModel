sirConcern := KEModelPart new.
sirConcern attributes: 
    {#status->#(#S #I #R)}.
sirConcern addParameters: {#lambda. #gamma}.
sirConcern
   addTransitionFrom: {(#status -> #S)}
   to: {(#status -> #I)}
   probability: 'lambda'.
sirConcern changeParameter: #lambda 
   value: 'taugen*I*alphagen'.
sirConcern
   addTransitionFrom: {(#status -> #I)} 
   to: {(#status -> #R)}
   probability: 'gamma'.
        
maConcern := KEModelPart new.
maConcern addParameter:  #alphagen
   value: 'deg'. 
maConcern addParameter:  #taugen
   value: 'beta'.    

model := sirConcern + maConcern.
    
model atCompartment: {(#status -> #S)} put: 0.975.
model atCompartment: {(#status -> #I)} put: 0.025.
model atParameter: #beta assignValue: 0.6.
model atParameter: #gamma assignValue: 1.
model atParameter: #deg assignValue: 5.

simulator := KESimulator
   new: #RungeKutta
   from: 0.0
   to: 15
   step: 0.1.
simulator executeOn: model.

chart := KEChart new.
chart addDataFrame: 
    (simulator timeSeriesOutputsAt: 
        {(#status -> #I)}).
chart  yLabel: 'Infectious'.
chart legendTitle: 'Mass Action'.
chart plot.
