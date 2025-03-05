//import kidney-logo from './kidney-logo.svg';
//import logo from '../public/kidneys_logo2.png';
import React from 'react';
import { Vitessce } from 'vitessce';
//import TimecoursePlot from './components/TimecoursePlot';
import { myViewConfig } from './my-view-config';

export default function App() {
  return (
    <Vitessce
      config={myViewConfig}
      theme="light"
    />
  );
}
