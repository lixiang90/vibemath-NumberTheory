import { renderToString } from "katex";
import { escape, href, inline, repository, sourceUrl } from "./content.mjs";

const library = "数论研究文库";
const home = "index.html";
const guide = "guide/index.html";

function shell(output, title, description, current, content) {
  const url = target => href(output, target);
  const nav = (target, label, section) => `<a href="${url(target)}"${current === section ? ' aria-current="page"' : ""}>${label}</a>`;
  return `<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="description" content="${escape(description)}">
  <title>${escape(title)} · VibeMath</title>
  <link rel="stylesheet" href="${url("assets/katex/katex.min.css")}">
  <link rel="stylesheet" href="${url("assets/site.css")}">
</head>
<body>
  <a class="skip" href="#main">跳转至正文</a>
  <header class="mast"><div class="wrap mast-inner">
    <a class="brand" href="${url(home)}" aria-label="VibeMath 数论研究文库首页" lang="en">Vibe<i>Math</i></a>
    <nav aria-label="主导航">
      ${nav(home, "研究目录", "home")}${nav(guide, "发表指南", "guide")}
      <a href="${repository}" lang="en">GitHub</a>
    </nav>
  </div></header>
  <main id="main" class="wrap">${content}</main>
  <footer class="wrap footer"><span>VibeMath <span class="footer-divider">/</span> 数论研究文库</span><a href="${repository}">源代码仓库</a></footer>
</body>
</html>
`;
}

function breadcrumbs(output, topic) {
  return `<nav class="breadcrumbs" aria-label="面包屑"><a href="${href(output, home)}">研究目录</a>${topic ? `<span aria-hidden="true">/</span><a href="${href(output, topic.output)}">${escape(topic.title)}</a>` : ""}</nav>`;
}

function featured(topics) {
  for (const topic of topics) {
    const paper = topic.papers.find(paper => paper.formula);
    if (!paper) continue;
    const formula = renderToString(paper.formula, { displayMode: true, throwOnError: true });
    return `<aside class="feature">
      <p class="feature-label">${paper.example ? "示例选读" : "论文选读"}<span> / </span>${escape(topic.title)}</p>
      <div class="feature-formula">${formula}</div>
      <a href="${href(home, paper.output)}">${inline(paper.title)}</a>
    </aside>`;
  }
  return "";
}

export function homeView(topics) {
  const rows = topics.map(topic => {
    const examples = topic.papers.filter(paper => paper.example).length;
    const count = topic.papers.length - examples;
    return `<li><a class="topic-row" href="${href(home, topic.output)}">
      <div><span class="directory" lang="en">${escape(topic.id)}</span><h3>${escape(topic.title)}</h3></div>
      <p>${escape(topic.description)}</p>
      <span class="counts">${count} 篇成果${examples ? `<small>${examples} 篇示例</small>` : ""}</span>
    </a></li>`;
  }).join("\n");
  const intro = "按研究领域收录数论论文与研究笔记，可在线阅读推导、公式与证明。";
  return shell(home, library, intro, "home", `
    <section class="home-intro"><div class="intro-copy"><h1>${library}</h1><p class="lead">${intro}</p></div>${featured(topics)}</section>
    <section aria-labelledby="topics"><div class="section-title"><h2 id="topics">研究目录</h2><span>${topics.length} 个领域</span></div>
      <ul class="topic-list">${rows}</ul>
    </section>`);
}

function entries(topic, papers, empty) {
  if (!papers.length) return `<p class="empty">${empty}</p>`;
  return `<ul class="entries">${papers.map(paper => {
    const credit = [paper.author, paper.date].filter(Boolean).map(escape).join(" · ");
    return `<li><div><h3><a href="${href(topic.output, paper.output)}">${inline(paper.title)}</a></h3>
      ${paper.abstract ? `<p>${inline(paper.abstract)}</p>` : ""}</div>
      ${credit ? `<span class="credit">${credit}</span>` : ""}</li>`;
  }).join("\n")}</ul>`;
}

export function topicView(topic) {
  const papers = topic.papers.filter(paper => !paper.example);
  const examples = topic.papers.filter(paper => paper.example);
  return shell(topic.output, topic.title, topic.description, "", `
    ${breadcrumbs(topic.output)}
    <header class="topic-intro"><p class="directory" lang="en">${escape(topic.id)}</p><h1>${escape(topic.title)}</h1><p class="lead">${escape(topic.description)}</p></header>
    <section class="collection"><div class="section-title"><h2>研究成果</h2><span>${papers.length} 篇</span></div>
      ${entries(topic, papers, "本领域尚未收录研究成果。")}
    </section>
    ${examples.length ? `<section class="collection"><div class="section-title"><h2>排版示例</h2><span>${examples.length} 篇</span></div>${entries(topic, examples)}</section>` : ""}`);
}

export function articleView(page, topic, { titleId, toc, body }) {
  const isGuide = !topic;
  const title = escape(page.title);
  const credit = [page.author, page.date].filter(Boolean).map(escape).join(" · ");
  const meta = isGuide ? "" : `<p class="article-kind">${page.example ? "排版示例" : "研究成果"}${credit ? `<span> / </span>${credit}` : ""}</p>`;
  const contents = toc.map(item => `<li><a href="#${encodeURIComponent(item.id)}">${escape(item.label)}</a></li>`).join("");
  return shell(page.output, page.title, page.abstract || page.title, isGuide ? "guide" : "", `
    ${breadcrumbs(page.output, topic)}
    <header class="article-heading">${meta}<h1 id="${escape(titleId)}">${inline(page.title)}</h1>
      ${page.abstract ? `<p class="abstract"><strong>摘要</strong>${inline(page.abstract)}</p>` : ""}
    </header>
    <div class="reading-layout">
      <aside class="contents">
        ${toc.length ? `<nav aria-label="本文目录"><h2>本文目录</h2><ul>${contents}</ul></nav>` : ""}
        <a class="source-link" href="${sourceUrl(page.source)}" aria-label="查看 Markdown 原文（${title}）">查看 Markdown 原文</a>
      </aside>
      <article class="prose">${body}</article>
    </div>`);
}
