/*
 * jQueryJmol - jQuery Plugin
 * Complete replacement for Jmol.js
 * 
 * Copyright (c) 2012 Gusts 'gusC' Kaksis
 * 
 * Version: 2.0.0 (10/08/2012)
 * Requires: jQuery v1.6+
 * Adapted to jQuery v2.0.0 by Miguel Vazquez (10/05/2013)
 *
 * Licensed under the MIT license:
 *   http://www.opensource.org/licenses/mit-license.php
 */
 
(function($){
	/*
	* Jmol Java applet's wrapper class
	* Used to expose external interface methods to JavaScript, when applet becomes ready
	*/
	var JmolWrapper = function(id, options){
		this.id = id;
		this.options = options;
		this._applet = null;
		this._callbackCounter = {};
		this._callbackNames = {
			animate : 'AnimFrameCallback',
			echo : 'EchoCallback',
			hover : 'HoverCallback',
			load : 'LoadStructCallback',
			measure : 'MeasureCallback',
			message : 'MessageCallback',
			minimize : 'MinimizationCallback',
			pick : 'PickCallback',
			resize : 'ResizeCallback',
			sync : 'SyncCallback', // callback not event!
			script : 'ScriptCallback'
		};
		this._scriptCache = [];
		
		// Wrapped methods
		
		/**
		* Get a property from Jmol
		* @param string - property name
		* @param string - property (sub)value
		* @return object - key:value pairs
		*/
		this.getProperty = function(property, value){
			var json = '[]';		
			if (this._applet){
				if (typeof value != 'undefined'){
					json = this._applet.getPropertyAsJSON(property, value);
				} else {
					json = this._applet.getPropertyAsJSON(property);
				}
			}
			return jQuery.parseJSON(json);
		};
		/**
		* Send a script to Jmol
		* @param string - script source
		* @return void
		*/
		this.script = function(source){
			if (this._applet){
				this._applet.script(source);
			} else {
				this._scriptCache.push(source);
			}
		};
		/**
		* Send a script to Jmol (bypases the queue and returns callback information)
		* @param string - script source
		* @return string
		*/
		this.scriptWait = function(source, statusParams){
			if (this._applet){
				if (typeof statusParams != 'undefined'){
					return this._applet.scriptWait(source, statusParams);
				} else {
					return this._applet.scriptWait(source);
				}
			}
			return '';
		};
		/**
		* Send a script to Jmol (bypases the queue and wait for any output from messages, echo, print, etc.)
		* @param string - script source
		* @return string
		*/
		this.scriptWaitOutput = function(source){
			if (this._applet){
				return this._applet.scriptWaitOutput(source);
			} 
			return '';
		};
		/**
		* Send a script to Jmol (put it on the queue)
		* @param string - script source
		* @return string
		*/
		this.scriptNoWait = function(source){
			if (this._applet){
				return this._applet.scriptNoWait(source);
			}
			return '';
		};
		/**
		* Check script syntax and report back with a string message
		* @param string - script source
		* @return string
		*/
		this.scriptCheck = function(source){
			if (this._applet){
				return this._applet.scriptCheck(source);
			}
			return '';
		};
		
		// Private methods
		
		/**
		* When applet calls a ready callback, this method must be invoked to store applet's 
		* interface internally. Also it will perform any cached scripts.
		* @param java - Jmol applet's interface
		*/
		this._ready = function(applet){
			this._applet = applet;
			if (this._scriptCache.length > 0){
				for (var a in this._scriptCache){
					this._applet.script(this._scriptCache[a]);
				}
			}
			this._scriptCache = [];
		};
		/**
		* Clean-up after destruction of Jmol
		*/
		this._destroy = function(){
			delete this._scriptCache;
			delete this._applet;
			delete this.options;
			delete this.id;
		};
		
		// Internal callback manager
		
		/**
		* Check and set up a callback if neccessary
		* @param string - event name without jmol_ part
		*/
		this._setCallback = function(type){
			if (typeof this._callbackCounter[type] == 'undefined'){
				this._callbackCounter[type] = 1;
				this.script('set ' + this._callbackNames[type] + ' "jQuery.jmol.' + type + '"');
			}
			this._callbackCounter[type] ++;
		};
		/**
		* Signal the manager about a removal of an event. If all the listeners are removed
		* also remove a callback from Jmol
		* @param string - event name without jmol_ part
		*/
		this._removeCallback = function(type){
			if (typeof this._callbackCounter[type] != 'undefined'){
				this._callbackCounter[type] --;
				if (this._callbackCounter[type] <= 0){
					this.script('set ' + this._callbackNames[type] + ' NONE');
				}
			}
		};
	};
	/**
	* Jmol jQuery plugin class
	*/
	var jQueryJmolPlugin = (function(){
		/**
		* Default option set
		*/
		var _defaults = {
			// Jmol initialization properties
			appletUrl : '', // URL of a directory where applet's jar file resides
			useSigned : false, // Use self signed version
			memLimit: 512, // Java memory limit in Megabytes
			width: 400, // Applets width in pixels
			height: 300, // Applets height in pixels
			menuUrl : '', // URL of a menu file
			modelUrl : '', // URL of an initial model file
			background: '#000000', // Background color
			
			// jQuery-Jmol callbacks
			onReady: null, // param: array of arguments
			onDestroy: null, // param: array of arguments
			onSync: null, // param: array of arguments
			onEval: null // param: array of arguments
		},
		/**
		* Unsigned applet file
		*/
		_appletFile = '/js-find/jmol/JmolApplet0.jar', 
		/**
		* Signed applet file
		*/
		_appletFileSigned = '/js-find/jmol/JmolAppletSigned0.jar', 
		/**
		* HTML template for Jmol applet
		*/
		_htmlTemplate = '<object type="application/x-java-applet" id="%id%" name="%name%" width="%width%" height="%height%"%add_attr%>'
			+ '<param name="syncId" value="%sync_id%"/>'
			+ '<param name="progressbar" value="true">'
			+ '<param name="progresscolor" value="blue">'
			+ '<param name="boxbgcolor" value="%bg_color%"/>'
			+ '<param name="boxfgcolor" value="black">'
			+ '<param name="boxmessage" value="Downloading JmolApplet ...">'
			+ '<param name="mayscript" value="mayscript">'
			+ '<param name="codebase" value="%applet_url%" />'
			+ '<param name="archive" value="%applet_file%" />'
			+ '<param name="code" value="JmolApplet.class" />'
			+ '<param name="java_arguments" value="%java_args%"/>'
			//+ '<param name="script" value="%script%"/>'

			+ '%add_param%'

			+ '<p>You do not have Java applets enabled in your web browser, or your browser is blocking this applet.<br>'
			+ 'Check the warning message from your browser and/or enable Java applets in<br>'
			+ 'your web browser preferences, or install the Java Runtime Environment from <a href="http://www.java.com">www.java.com</a><br></p>'
			+ '</object>',
		/**
		* We set thease callbacks only if options are set
		* so we don't bring unnecessary load to JavaScript
		* 
		* Also thease callbacks are set after Jmol has initialized
		*/
		
		/**
		* Thease callbacks have to be set before the applet has initialized (documentation says so :)
		*/
		_cbBefore = {
			eval : '<param name="evalCallback" value="$.jmol.eval" />',
			ready : '<param name="appletReadyCallback" value="$.jmol.ready" />'
		},
		/**
		* Java applet Class ID
		*/
		_windowsClassId = 'clsid:8AD9C840-044E-11D1-B3E9-00805F499D93',
		/**
		* Windows CAB URL for Java installer
		*/
		_windowsCabUrl = 'http://java.sun.com/update/1.6.0/jinstall-6u22-windows-i586.cab',
		/**
		* Internal applet counter, for unique ID generation
		*/
		appletCounter = 0;
		
		/**
		* Main entry point for jQuery plugin initialization
		* @param mixed - object for initialization options, string internal commands (hide, show, destroy)
		* @return jQuery
		*/
		var _process = function(command){
			return this.each(function(i, item) {
				var $item = $(item);
				if ($item.data('jmol')){
					var jmol = $item.data('jmol');
					if (typeof command == 'string'){
						// Perform some jQuery or DOM related tasks with Jmol Applet
						switch (command){
							case 'hide':
								// Hide Jmol Applet
								// We can't use dislpay:none, it will break Java
								$item.find('object').css('width', '2px');
								$item.find('object').css('height', '2px');
								break;
							case 'show':
								// Restore Jmol Applet in to the view
								$item.find('object').css('width', jmol.options['width'] + 'px');
								$item.find('object').css('height', jmol.options['height'] + 'px');
								break;
							case 'destroy':
								// Remove an applet from the view completely and forget it
								jmol._destroy();
								$item.find('object').remove();
								$item.data('jmol', null);
								break;
							case 'option':
								// Update single or multiple options
								var prop = arguments[1];
								if (typeof prop == 'string'){
									var val = arguments[2];
									if (typeof val != 'undefined'){
										_updateOption($item, jmol, prop, val);
									}
								} else if (typeof prop == 'object'){
									_updateOptionObj($item, jmol, prop);
								}
								break;
						}
					} else if (typeof command == 'object'){
						// Update options
						_updateOptionObj($item, jmol, command);
					}
				} else {
					// This will allow Jmol to initialize with default options
					var options = $.extend({}, _defaults, {});
					if (typeof command == 'object'){
						options = $.extend(options, command || {});
					}
					appletCounter ++;
					var id = $item.attr('id');
					if (typeof id == 'undefined'){
						// We need an unique ID, so here we generate it if there is none
						id = 'jmolApplet' + appletCounter
						$item.attr('id', id);
					}
					var jmol = new JmolWrapper(id, options);
					// Load default assets
					if (options['menuUrl'] != null && options['menuUrl'].length > 0){
						jmol.script('load MENU ' + options['menuUrl']);
					}
					if (options['modelUrl'] != null && options['modelUrl'].length > 0){
						jmol.script('load ' + options['modelUrl']);
					}
					// Set callback functions
					if (options['onSync'] != null){
						jmol._setCallback('sync');
					}
					// Insert HTML block
					$item.html(_appletBuildHtml('_' + id, options));
					// Mark Jmol initialized
					$item.data('jmol', jmol);
				}
			});
		},
		/**
		* Update option array
		*/
		_updateOptionObj = function($item, jmol, options){
			for (var a in options){
				_updateOption($item, jmol, a, options[a]);
			}
		},
		/**
		* Update single option
		*/
		_updateOption = function($item, jmol, name, value){
		  switch (name){
				case 'background':
					jmol.script('background ' + value.replace('#', 'x'));
					break;
				case 'menuUrl':
					jmol.script('load MENU "' + value + '"');
					break;
				case 'modelUrl':
					jmol.script('load "' + value + '"');
					break;
				case 'onSync':
					jmol.setCallback('sync');
					break;
				case 'width':
					$item.find('object').css('width', value + 'px');
					break;
				case 'height':
					$item.find('object').css('width', value + 'px');
					break;	
			}
			jmol.options[name] = value;
		},	
		/**
		* Build applet's HTML block from a template
		* @param string - ID attribute
		* @param string - additional class attributes
		* @param object - jquery plugin options
		* @return string - applets HTML code
		*/
		_appletBuildHtml = function(id, options){
			var add_attr = '';
			if (navigator.userAgent){
				if (navigator.userAgent.indexOf('MSIE') != -1){
					// IE - add classid and codebase
					add_attr = ' classid="' + _windowsClassId + '" codebase="' + _windowsCabUrl + '"';
				}
			}
			var add_param = '';
			if (navigator.platform){
				if (navigator.platform.indexOf('Mac') != -1){
					// MacOS - add command thread to overcome some Java security restrictions
					// like java.security.AccessControlException: access denied (java.net.SocketPermission...)
					add_param = '<param name="UseCommandThread" value="true">'; 
				}
			}
			// We need ready internaly, so this goes on by default
			add_param += _cbBefore['ready'];
			// Eval is optional, but is required before initialization
			if (options['onEval'] !== null){
				add_param += _cbBefore['eval'];
			}
			
			var html = _htmlTemplate.replace('%add_attr%', add_attr);
			html = html.replace('%add_param%', add_param);
			html = html.replace('%sync_id%', ("" + Math.random()).substring(3));
			html = html.replace('%id%', id);
			html = html.replace('%name%', id);
			html = html.replace('%width%', options['width']);
			html = html.replace('%height%', options['height']);
			html = html.replace('%applet_url%', options['appletUrl']);
			html = html.replace('%applet_file%', (options['useSigned'] ? _appletFileSigned : _appletFile));
			html = html.replace('%java_args%', '-Xmx' + options['memLimit'] + 'm');
			html = html.replace('%bg_color%', options['background']);
			return html;
		},
		/**
		* Process calllbacks coming from Jmol
		* @param string - callback name
		* @param object - list of arguments passed to listener
		* @return integer - for sync callback
		*/
		_callback = function(type, args){
			var shift = Array.prototype.shift;
			var unshift = Array.prototype.unshift;
			var id = shift.call(args);
			var $item = $('#' + id.substr(1));
			if ($item.data('jmol')){
				var jmol = $item.data('jmol');
				switch (type){
					case 'ready':
						if (args[1]){
							jmol._ready(args[2]);
							if (jmol.options['onReady'] != null){
								jmol.options.onReady(jmol);
							}
						} else {
							if (jmol.options['onDestroy'] != null){
								jmol.options.onDestroy();
							}
						}
						break;
					case 'eval':
						if (jmol.options['onEval'] != null){
							jmol.options.onEval(jmol, args);
						}
						break;
					case 'sync':
						if (jmol.options['onSync'] != null){
							return jmol.options.onSync(jmol, args);
						}
						return 1
						break;
				}
			}
		},
		/**
		* Process events coming from Jmol
		* @param string - callback name
		* @param object - list of arguments passed to listener
		*/
		_event = function(type, args){
			var shift = Array.prototype.shift;
			var unshift = Array.prototype.unshift;
			var id = shift.call(args);
			var $item = $('#' + id.substr(1));
			if ($item.data('jmol')){
				var jmol = $item.data('jmol');
				var e = $.Event('jmol_' + type);
				var eargs = [jmol];
				if (type == 'animate'){
					eargs.push({
						frameIdx : args[0],
						fileNum : args[1],
						modelNum : args[2],
						firstFrame : args[3],
						lastFrame : args[4],
						animState : args[5],
						animDir : args[6],
						direction : args[7]
					});
				} else if (type == 'hover' || type == 'pick'){
					eargs.push({
						label : args[0],
						atomIdx : args[1]
					});
				} else if (type == 'load'){
					eargs.push({
						fileUrl : args[0],
						fileName : args[1],
						title : args[2],
						errorMsg : args[3],
						status : args[4],
						frameBefore : args[5],
						frameLast : args[6]
					});
				} else if (type == 'minimize'){
					eargs.push({
						status : args[0],
						iteration : args[1],
						energy : args[2],
						energyDelta : args[3]
					});
				} else if (type == 'measure'){
					eargs.push({
						label : args[0],
						status : args[2],
						value : args[3]
					});
				} else if (type == 'echo' || type == 'message' || type == 'script'){
					eargs.push(args[0]);
				} else if (type == 'resize'){
					eargs.push({
						width : args[0],
						height : args[0]
					});
				}
				$item.triggerHandler(e, eargs);
			}
		};
		
	//	// We need to override on() and off() methods
	//	var _oldBind = $.event.add;
	//	$.event.add = function(elem, types, handler, data){
	//		$(elem).each(function(i, item){
	//			if ($(item).data('jmol')){
	//				var jmol = $(item).data('jmol');
	//				var etypes = types.split(' ');
	//				for (var a in etypes){
	//					switch (etypes[a]){
	//						case 'jmol_animate':
	//						case 'jmol_echo':
	//						case 'jmol_hover':
	//						case 'jmol_load':
	//						case 'jmol_measure':
	//						case 'jmol_message':
	//						case 'jmol_minimize':
	//						case 'jmol_pick':
	//						case 'jmol_resize':
	//						case 'jmol_script':
	//							jmol._setCallback(etypes[a].substr(5));
	//				  		break;
	//					}
	//				}
	//			}
	//		});
	//		// Pass the controll back to old one
	//		return _oldBind(elem, types, handler, data);
	//	}

	//	var _oldUnBind = $.event.remove;

  //// Changed out by Miguel Vazquez
  ////   I don't know what this does, but it breaks my jquery. I changed fn for
  ////   event for simetry with the previous. But it seems to have no effect even
  ////   to comment it all out.
  ////
	//	//$.fn.remove = function(elem, types, handler, pos){
	//	$.event.remove = function(elem, types, handler, pos){
  // 
	//		$(elem).each(function(i, item){
	//			if ($(item).data('jmol')){
	//				var jmol = $(item).data('jmol');
	//				var etypes = types.split(' ');
	//				for (var a in etypes){
	//					switch (etypes[a]){
	//						case 'jmol_animate':
	//						case 'jmol_echo':
	//						case 'jmol_hover':
	//						case 'jmol_load':
	//						case 'jmol_measure':
	//						case 'jmol_message':
	//						case 'jmol_minimize':
	//						case 'jmol_pick':
	//						case 'jmol_resize':
	//						case 'jmol_script':
	//				  		jmol._removeCallback(etypes[a].substr(5));
	//				  		break;
	//					}
	//				}
	//			}
	//		});
	//		// Pass the controll back to old one
	//		return _oldUnBind(elem, types, handler, pos);
	//	}
		
		// Export methods to public space
		return {
			process : _process,
			callbacks : {
				eval : function(){
					_callback('eval', arguments);
				},
				sync : function(){
					return _callback('sync', arguments);
				},
				message : function(){
					_event('message', arguments);
				},
				echo : function(){
					_event('echo', arguments);
				},
				script : function(){
					_event('script', arguments);
				},
				ready : function(){
					_callback('ready', arguments);
				},
				load : function(){
					_event('load', arguments);
				},
				hover : function(){
					_event('hover', arguments);
				},
				pick : function(){
					_event('pick', arguments);
				},
				measure : function(){
					_event('measure', arguments);
				},
				animate : function(){
					_event('animate', arguments);
				},
				minimize : function(){
					_event('minimize', arguments);
				},
				resize : function(){
					_event('resize', arguments);
				}
			}
		};
	})();
	
	// Register plugin
	$.fn.extend({
		// register jQuery-Jmol plugin
		jmol : jQueryJmolPlugin.process
	});
	
	// Register callback methods into global jQuery namespace
	$.extend({
		jmol : jQueryJmolPlugin.callbacks
	});
	
})($);
