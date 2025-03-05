import { registerPluginView } from '@vitessce/plugin-api';
import TimecoursePlot from './components/TimecoursePlot';

export function registerVitesscePlugins() {
    registerPluginView('TimecoursePlotView', TimecoursePlot);
}
