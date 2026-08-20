/* open-king documentation: theme toggle, copy buttons, heading anchors,
   table of contents. Vanilla JS, no dependencies, no network access. */

(function () {
  "use strict";

  var STORAGE_KEY = "open-king-theme";
  var root = document.documentElement;

  /* ---------------------------------------------------------------- theme */

  function systemTheme() {
    return window.matchMedia &&
      window.matchMedia("(prefers-color-scheme: dark)").matches
      ? "dark"
      : "light";
  }

  function activeTheme() {
    var explicit = root.getAttribute("data-theme");
    return explicit === "light" || explicit === "dark" ? explicit : systemTheme();
  }

  function describeToggle(button) {
    var next = activeTheme() === "dark" ? "light" : "dark";
    var label = "Switch to " + next + " theme";
    button.setAttribute("aria-label", label);
    button.setAttribute("title", label);
  }

  function initTheme() {
    var button = document.querySelector(".theme-toggle");
    if (!button) {
      return;
    }
    describeToggle(button);

    button.addEventListener("click", function () {
      var next = activeTheme() === "dark" ? "light" : "dark";
      root.setAttribute("data-theme", next);
      try {
        localStorage.setItem(STORAGE_KEY, next);
      } catch (err) {
        /* storage unavailable: the choice lasts for this page only */
      }
      describeToggle(button);
    });

    if (window.matchMedia) {
      var query = window.matchMedia("(prefers-color-scheme: dark)");
      var onChange = function () {
        describeToggle(button);
      };
      if (query.addEventListener) {
        query.addEventListener("change", onChange);
      } else if (query.addListener) {
        query.addListener(onChange);
      }
    }
  }

  /* ------------------------------------------------------------ narrow nav */

  function initNav() {
    var button = document.querySelector(".nav-toggle");
    var sidebar = document.querySelector(".sidebar");
    if (!button || !sidebar) {
      return;
    }

    button.addEventListener("click", function () {
      var open = sidebar.getAttribute("data-nav-open") === "true";
      sidebar.setAttribute("data-nav-open", open ? "false" : "true");
      button.setAttribute("aria-expanded", open ? "false" : "true");
    });

    document.addEventListener("keydown", function (event) {
      if (event.key === "Escape") {
        sidebar.setAttribute("data-nav-open", "false");
        button.setAttribute("aria-expanded", "false");
      }
    });
  }

  /* --------------------------------------------------------- copy buttons */

  function copyText(text) {
    if (navigator.clipboard && navigator.clipboard.writeText) {
      return navigator.clipboard.writeText(text);
    }
    return new Promise(function (resolve, reject) {
      var field = document.createElement("textarea");
      field.value = text;
      field.setAttribute("readonly", "");
      field.style.position = "absolute";
      field.style.left = "-9999px";
      document.body.appendChild(field);
      field.select();
      var ok = false;
      try {
        ok = document.execCommand("copy");
      } catch (err) {
        ok = false;
      }
      document.body.removeChild(field);
      if (ok) {
        resolve();
      } else {
        reject(new Error("copy failed"));
      }
    });
  }

  function initCopyButtons() {
    var blocks = document.querySelectorAll(".code-block");
    Array.prototype.forEach.call(blocks, function (block) {
      if (block.classList.contains("code-block--plain")) {
        return;
      }
      var code = block.querySelector("pre code");
      if (!code) {
        return;
      }

      var button = document.createElement("button");
      button.type = "button";
      button.className = "copy-button";
      button.textContent = "copy";
      button.setAttribute("aria-label", "Copy code to clipboard");

      button.addEventListener("click", function () {
        copyText(code.textContent).then(
          function () {
            button.textContent = "copied";
            button.setAttribute("data-copied", "true");
            window.setTimeout(function () {
              button.textContent = "copy";
              button.removeAttribute("data-copied");
            }, 1600);
          },
          function () {
            button.textContent = "select and copy";
            window.setTimeout(function () {
              button.textContent = "copy";
            }, 2400);
          }
        );
      });

      block.appendChild(button);
    });
  }

  /* ------------------------------------------------------ heading anchors */

  function slug(text) {
    return text
      .toLowerCase()
      .replace(/[^\w\s.-]/g, "")
      .trim()
      .replace(/[\s.]+/g, "-")
      .replace(/-+/g, "-");
  }

  function headingLabel(heading) {
    var copy = heading.cloneNode(true);
    var anchor = copy.querySelector(".heading-anchor");
    if (anchor) {
      copy.removeChild(anchor);
    }
    return copy.textContent.replace(/\s+/g, " ").trim();
  }

  function initHeadings(headings) {
    var used = {};
    Array.prototype.forEach.call(headings, function (heading) {
      if (!heading.id) {
        var base = slug(heading.textContent) || "section";
        var id = base;
        var n = 2;
        while (used[id] || document.getElementById(id)) {
          id = base + "-" + n;
          n += 1;
        }
        heading.id = id;
      }
      used[heading.id] = true;

      if (heading.querySelector(".heading-anchor")) {
        return;
      }
      var anchor = document.createElement("a");
      anchor.className = "heading-anchor";
      anchor.href = "#" + heading.id;
      anchor.textContent = "#";
      anchor.setAttribute("aria-label", "Permalink to this section");
      heading.appendChild(anchor);
    });
  }

  /* -------------------------------------------------------------- contents */

  function initToc(headings) {
    var toc = document.getElementById("toc");
    if (!toc || headings.length < 3) {
      return;
    }

    var title = document.createElement("p");
    title.className = "toc-title";
    title.textContent = "On this page";

    var list = document.createElement("ul");
    list.className = "toc-list";

    var links = [];
    Array.prototype.forEach.call(headings, function (heading) {
      var item = document.createElement("li");
      item.className =
        heading.tagName === "H3" ? "toc-item toc-item--sub" : "toc-item";

      var link = document.createElement("a");
      link.className = "toc-link";
      link.href = "#" + heading.id;
      link.textContent = headingLabel(heading);

      item.appendChild(link);
      list.appendChild(item);
      links.push(link);
    });

    toc.appendChild(title);
    toc.appendChild(list);
    toc.hidden = false;

    if (!window.IntersectionObserver) {
      return;
    }

    var visible = {};
    var observer = new window.IntersectionObserver(
      function (entries) {
        entries.forEach(function (entry) {
          visible[entry.target.id] = entry.isIntersecting;
        });

        var current = null;
        Array.prototype.forEach.call(headings, function (heading) {
          if (visible[heading.id] && !current) {
            current = heading.id;
          }
        });

        links.forEach(function (link) {
          var match = current && link.hash === "#" + current;
          link.classList.toggle("toc-link--current", Boolean(match));
        });
      },
      { rootMargin: "0px 0px -70% 0px", threshold: 0 }
    );

    Array.prototype.forEach.call(headings, function (heading) {
      observer.observe(heading);
    });
  }

  /* ------------------------------------------------------------------ boot */

  function boot() {
    var article = document.querySelector(".prose");
    var headings = article ? article.querySelectorAll("h2, h3") : [];

    initTheme();
    initNav();
    initCopyButtons();
    initHeadings(headings);
    initToc(headings);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", boot);
  } else {
    boot();
  }
})();
